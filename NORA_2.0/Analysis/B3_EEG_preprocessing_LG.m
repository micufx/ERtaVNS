clc; clear; close all;

%%
% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% EEG preprocessing (Local-Global | corrected trigger workflow)
% This code uses the ICA learned previously (from concatenated E1+E2),
% performs ERP-prep filtering, applies the corrected Local-Global labels,
% epochs around the fifth-tone onset, baselines, applies session-level
% ICLabel auto-rejection, then per-experiment epoch rejection, interpolation,
% re-referencing, plotting, and saving.
%
% CURRENT 4-bit acquisition trigger architecture:
%   S5..S8 = T1..T4; S9 = T5 anchor.
%   Immediately after S9: S1=FLS, S2=FLD, S3=RLS, S4=RLD.
%   On paired-stimulation trials S10 follows the condition marker.
%   S10 labels selected GD and matched GS trials; it is not the epoch anchor.
%
% label_LG_events.m uses S9 as the zero-time T5 event, reads the following
% S1..S4 condition marker and optional S10, then annotates the S9 event in place.
% No marker latency is shifted. Habituation is identified from the first 25
% trials of each 205-trial block and remains available in lg_phase.
%
% set_LG_epoch_labels.m then assigns one scalar trial-defining condition to
% each epoch from the T5 anchor at 0 ms, avoiding ambiguity from neighboring
% events that also fall inside the long epoch window.
%
% IMPORTANT:
%   The corrected acquisition provides an onset-aligned S9 T5 anchor followed
%   by an explicit S1..S4 condition label. No post-hoc latency shifting or
%   block-based condition reconstruction is used.
%
% Miguel Contreras-Altamirano, updated 2026

%% Paths & files

num_subject = 1;    % Start analyzing from participant:

analysisVisit = 1;   % <-- change to 2 when processing Visit 2
PROJECT_ROOT = fileparts(LG_THIS_SCRIPT_DIR);
mainpath = 'L:\Downloads\eeglab2023.0\'; % EEGLAB folder (edit if needed)
path     = fullfile(PROJECT_ROOT,'Experiments','Data');
outpath  = fullfile(PROJECT_ROOT,'Analysis_Output',sprintf('Visit%d',analysisVisit));

% Add EEGLAB properly
if exist('eeglab','file') ~= 2
    addpath(mainpath);
end

% --- Find participants that have a concatenated ICA set: out\<Sxx>\ica_<Sxx>.set ---
dSubs = dir(outpath);
isSub = [dSubs.isdir] & ~ismember({dSubs.name},{'.','..'});
subNames = {dSubs(isSub).name};
subNames = subNames(~cellfun('isempty', regexp(subNames, '^S\d+', 'once')));         % keep S##

% Build list of participants with the expected concatenated ICA file
participant_list = {};
for si = 1:length(subNames)
    p  = subNames{si};
    fC = fullfile(outpath, p, ['ica_' p '.set']);
    if exist(fC,'file')
        participant_list{end+1} = p; %#ok<SAGROW>
    end
end
if isempty(participant_list)
    error('No concatenated ICA sets found (expected out\\Sxx\\ica_Sxx.set).');
end
participant_list = sort(participant_list);

%% Parameters

% ERP parameters (aligned with classical ERP prep; Bekinschtein used 0.5–20 Hz)
HPF_erp     = 0.5;               % High-pass (Hz)
LPF_erp     = 20;                % Low-pass  (Hz)

epochWin    = [-1.4 0.7];          % s, time-locked to 5th-tone trial onset marker
% Baseline is computed dynamically as -800..-600 ms relative to T5.

% Four elementary LG conditions. Epoching directly on these labels keeps the
% full 2 x 2 Local x Global structure available throughout preprocessing.
condLabels = {'FLS','RLS','FLD','RLD'};

% Session-level ICLabel thresholds. Any component with at least %num confidence of being an artifact or
% less than %num confidence for brain
thr.muscle   = 0.30; % Muscle
thr.eye      = 0.90; % Eye
thr.heart    = 0.90; % Heart
thr.line     = 0.90; % Line Noise
thr.channel  = 0.50; % Channel Noise
thr.other    = 0.50; % Other
thr.brainMin = NaN; % Brain

% Per-experiment automatic epoch rejection (peak-to-peak μV); set [] to skip
autoRejPP = 120;

% Optional re-reference
re_reference = false;    % true → TP9/TP10; false → average (consistent with many ERP preps)

% Four selectable channels: each row produces an old AND a new figure.
% First column = Local (left in new / top in old); second = Global (right / bottom).
channelPairs = {'FZ','PZ'; 'CZ','C4'};
toneOnsets_ms = [-600 -450 -300 -150 0]; % T1..T5; task traces only
auxChannel = 'StimVNS'; % Confirmed recorded ear-electrode label
auxChannelAliases = {'StimtaVNS'}; % Optional explicit alias used in earlier configuration
auxPreprocessing = 'erp'; % 'erp' = ERP filters + baseline; 'raw' = neither (diagnostic)
auxRequired = true; % Stop on missing/ambiguous AUX instead of silently omitting it
auxLineWidth = 4.8; % Bold red and lilac auxiliary mean curves
fprintf('Running AUX FIX 2026-09-29: %s\n',mfilename('fullpath'));
assert(ismember(lower(string(auxPreprocessing)),["erp","raw"]),'auxPreprocessing must be erp or raw.');
auxUnitScaleToMicrovolts = 1; % Assumes imported AUX units are uV; change ONLY for a known unit conversion
auxControlColor = [0.85 0.05 0.12]; % red: retained control-run trials (S10-locked)
auxPairedColor = [0.65 0.30 0.85]; % lilac: ALL retained paired task trials (T5-locked)
expandPlotLimitsForAux = true; % Expand BOTH new panels equally to show AUX in actual uV
% AUX is captured before channel rejection/ICA/scalp rereferencing, then ERP
% filtered/baselined and matched to exactly the retained scalp epochs.
assert(iscell(channelPairs) && isequal(size(channelPairs),[2 2]),'Specify two channel pairs (2-by-2 cell array).');
assert(numel(toneOnsets_ms)==5 && toneOnsets_ms(5)==0 && all(diff(toneOnsets_ms)>0),'Specify ordered T1..T5 times ending at zero.');

% New control/figure settings. Existing task preprocessing remains unchanged.
controlBaseWin_ms = [-800 -600]; % Same usual task baseline, relative to S10 onset.
erpPlotLimits = [-8 8];          % Shared Local/Global ERP limits (microvolts).
topoPlotLimits = [-2 2];        % Shared colorbar limits (microvolts); [] = shared autoscale.
topoLocalWin_ms = [120 220];    % ALL Local maps, including Pure taVNS / Pure Sham.
topoGlobalWin_ms = [300 600];   % ALL Global maps, including Pure taVNS / Pure Sham.
topoInsetTopGap = 0.008;        % Gap below panel top, as a fraction of panel height.
processControl = true;

% Optional movie

%% Summary holders (LG; session-level + per-experiment)

Bad_components_concatenated = containers.Map('KeyType','char','ValueType','any');
IC_labels_concatenated      = containers.Map('KeyType','char','ValueType','char');

Bad_trials_E1        = containers.Map('KeyType','char','ValueType','any');
Num_bad_trials_E1    = containers.Map('KeyType','char','ValueType','double');
Num_good_trials_E1   = containers.Map('KeyType','char','ValueType','double');
nStd_E1              = containers.Map('KeyType','char','ValueType','double');
nDev_E1              = containers.Map('KeyType','char','ValueType','double');

Bad_trials_E2        = containers.Map('KeyType','char','ValueType','any');
Num_bad_trials_E2    = containers.Map('KeyType','char','ValueType','double');
Num_good_trials_E2   = containers.Map('KeyType','char','ValueType','double');
nStd_E2              = containers.Map('KeyType','char','ValueType','double');
nDev_E2              = containers.Map('KeyType','char','ValueType','double');

%% EEGLAB + Script 1 table

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% Loop participants

for sub = num_subject:length(participant_list)

    load(fullfile(outpath, 'Info_EEG_LG.mat'), 'T_LG');   % Load Script 1 table (bad channels per participant)

    % ===== Ensure output columns exist in T_LG (do this once, before looping) =====
    if ~ismember('Bad_components_concatenated', T_LG.Properties.VariableNames)
        T_LG.Bad_components_concatenated = repmat({[]}, height(T_LG), 1);
    end
    if ~ismember('IC_labels_concatenated', T_LG.Properties.VariableNames)
        T_LG.IC_labels_concatenated = repmat({''}, height(T_LG), 1);
    end

    if ~ismember('Bad_trials_E1', T_LG.Properties.VariableNames)
        T_LG.Bad_trials_E1 = repmat({[]}, height(T_LG), 1);
    end
    if ~ismember('Num_bad_trials_E1', T_LG.Properties.VariableNames)
        T_LG.Num_bad_trials_E1 = nan(height(T_LG), 1);
    end
    if ~ismember('Num_good_trials_E1', T_LG.Properties.VariableNames)
        T_LG.Num_good_trials_E1 = nan(height(T_LG), 1);
    end
    if ~ismember('nStd_E1', T_LG.Properties.VariableNames)
        T_LG.nStd_E1 = nan(height(T_LG), 1);
    end
    if ~ismember('nDev_E1', T_LG.Properties.VariableNames)
        T_LG.nDev_E1 = nan(height(T_LG), 1);
    end

    if ~ismember('Bad_trials_E2', T_LG.Properties.VariableNames)
        T_LG.Bad_trials_E2 = repmat({[]}, height(T_LG), 1);
    end
    if ~ismember('Num_bad_trials_E2', T_LG.Properties.VariableNames)
        T_LG.Num_bad_trials_E2 = nan(height(T_LG), 1);
    end
    if ~ismember('Num_good_trials_E2', T_LG.Properties.VariableNames)
        T_LG.Num_good_trials_E2 = nan(height(T_LG), 1);
    end
    if ~ismember('nStd_E2', T_LG.Properties.VariableNames)
        T_LG.nStd_E2 = nan(height(T_LG), 1);
    end
    if ~ismember('nDev_E2', T_LG.Properties.VariableNames)
        T_LG.nDev_E2 = nan(height(T_LG), 1);
    end

    participant = participant_list{sub};  % e.g., 'S01'

    % ===== Row in T_LG for THIS participant (must be inside the loop) =====
    rowIdx = find(strcmpi(T_LG.Subject_ID, participant), 1);
    if isempty(rowIdx)
        warning('No T_LG row found for %s; will skip table writes.', participant);
        badchans = [];
    else
        badchans = T_LG.Bad_channels{rowIdx};
        if isempty(badchans), badchans = []; end
    end

    out_subfold = fullfile(outpath, participant);
    inConcat    = fullfile(out_subfold, ['ica_' participant '.set']);

    fprintf('\n[%d/%d] %s\n', sub, length(participant_list), participant);

    %% --- Load concatenated ICA set (session-level; no epoching here) ---
    try
        EEGc = pop_loadset(inConcat);
    catch
        warning('pop_loadset failed for %s. Retrying once...', inConcat);
        pause(0.5);
        EEGc = pop_loadset(inConcat);
    end
    EEGc = eeg_checkset(EEGc);

    % (Optional) Inspect all ICs visually (topos)
    % pop_topoplot(EEGc, 0, 1:size(EEGc.icawinv,2), 'all-ica',[6 6], 0, 'electrodes','on');

    [EEGc,rej,lbltxt] = classify_eeg_components(EEGc,thr,out_subfold,participant);
    Bad_components_concatenated(participant) = rej;
    IC_labels_concatenated(participant) = char(strjoin(lbltxt,', '));

    % --- Write session-level results directly into T_LG (correct row) ---
    if ~isempty(rowIdx)
        T_LG.Bad_components_concatenated{rowIdx} = rej;
        T_LG.IC_labels_concatenated{rowIdx}      = char(strjoin(lbltxt, ', '));
    end

    % Save progress so stopping the script won't lose previous subjects
    save(fullfile(outpath,'Info_EEG_LG.mat'), 'T_LG');


    % Prepare control once, then reuse its retained epochs in E1/E2 figures.
    EEGcontrol=[];
    if processControl
        try
            EEGcontrol=preprocess_control_eeg(path,out_subfold,participant,analysisVisit, ...
                mainpath,HPF_erp,LPF_erp,epochWin,controlBaseWin_ms,thr,autoRejPP,re_reference,auxChannel,auxUnitScaleToMicrovolts,auxChannelAliases,auxPreprocessing,auxRequired);
        catch ME
            warning('B3:ControlFailed','%s V%d control unavailable: %s',participant,analysisVisit,ME.message);
        end
    end

    %% --- Process each experiment separately (E1 then E2) ---
    for eidx = 1:2

        % Loading experiments (E1 & E2)
        thisExp  = sprintf('E%d', eidx);                          % 'E1' or 'E2'
        subFolder=fullfile(path,participant,sprintf('Visit%d',analysisVisit));
        dVhdr=dir(fullfile(subFolder,sprintf('%sV%d%s*.vhdr',participant,analysisVisit,thisExp)));
        dVhdr=prefer_exact_stem(dVhdr,sprintf('%sV%d%s.vhdr',participant,analysisVisit,thisExp));
        if isempty(dVhdr)
            fprintf('   [%s V%d] No %sV%d%s*.vhdr. Skipping.\n',participant,analysisVisit,participant,analysisVisit,thisExp);
            continue;
        end

        %% Load raw experiment
        EEG = pop_loadbv(subFolder, dVhdr(1).name, [], []);
        EEG = eeg_checkset(EEG);

        %% Apply corrected trigger labels directly (NO latency shift)
        % S9 is the raw T5 anchor. label_LG_events reads the immediately
        % following S1..S4 condition label and optional S10, then relabels
        % the original S9 event in place without shifting its latency.
        [EEG, triggerReport] = label_LG_events(EEG, 'verbose', true);
        EEG = eeg_checkset(EEG, 'eventconsistency');

        [EEG, timing, T_LG] = check_effective_timing_LG( ...
            EEG, participant, T_LG, thisExp, ...
            'baselineBeforeX1_ms', 200, ...
            'postTone5_ms', 700, ...
            'expectedSOA_ms', 150, ...
            'expectedITI_ms', 2650, ...
            'store', 'median', ...
            'verbose', true);

        baseWin = timing.suggestedBaseWin_ms;

        %% Online reference note
        origRef = 'FCz';   % <- Used online
        EEG.ref = origRef; % 1) Tag in dataset
        EEG.etc.initial_reference.label  = origRef;   % 2) Save record in etc
        EEG.etc.initial_reference.kind   = 'hardware';
        EEG.etc.initial_reference.source = 'BrainVision online reference';
        EEG.comments = char(EEG.comments, ...
            sprintf('[%s] Initial hardware reference: %s', datestr(now,'yyyy-mm-dd HH:MM'), origRef));
        EEG = eeg_checkset(EEG);

        % Save a separate raw AUX stream; scalp processing below is unchanged.
        EEG=tag_aux_anchors(EEG);
        EEGauxRaw=extract_aux_channel(EEG,auxChannel,auxUnitScaleToMicrovolts,auxChannelAliases,auxRequired,fullfile(subFolder,dVhdr(1).name));
        taskAux=prepare_aux_epochs(EEGauxRaw,EEG,HPF_erp,LPF_erp,epochWin,baseWin,condLabels,auxPreprocessing);
        clear EEGauxRaw;

        %% Channel locations
        EEG = pop_chanedit(EEG, 'lookup', [mainpath 'plugins\dipfit5.4\standard_BESA\standard-10-5-cap385.elp']);
        EEG = eeg_checkset(EEG, 'eventconsistency');

        %% Bad channels removal. Remove the SAME channels as those removed before learning ICA (Script 1)
        if ~isempty(badchans)
            EEG = pop_select(EEG, 'rmchannel', badchans);
            EEG = eeg_checkset(EEG);
        end

        if EEG.nbchan ~= EEGc.nbchan % Safety check: channels must match the concatenated set used to learn ICA
            warning('   [%s %s] Channel count mismatch (this=%d, concat=%d). Skipping.', ...
                participant, thisExp, EEG.nbchan, EEGc.nbchan);
            continue;
        end

        %% Transfer ICA weights from concatenated set to this experiment
        EEG.icaact      = [];
        EEG.icawinv     = EEGc.icawinv;
        EEG.icasphere   = EEGc.icasphere;
        EEG.icaweights  = EEGc.icaweights;
        EEG.icachansind = EEGc.icachansind;
        EEG = eeg_checkset(EEG);

        %% ERP filtering (0.5-20 Hz), epoch on FLS/RLS/FLD/RLD & baseline
        EEG_erp = pop_eegfiltnew(EEG, 'locutoff', HPF_erp);  % High-pass
        EEG_erp = pop_eegfiltnew(EEG_erp, 'hicutoff', LPF_erp); % Low-pass
        EEG_erp = eeg_checkset(EEG_erp);

        EEG_erp = pop_epoch(EEG_erp, condLabels, epochWin, 'epochinfo','yes'); % T5 onset anchors
        EEG_erp = pop_rmbase(EEG_erp, baseWin);                                % Baseline correction
        EEG_erp = eeg_checkset(EEG_erp);

        % IMPORTANT: long epochs also contain markers from neighboring patterns.
        % Force the LG condition fields to refer ONLY to the T5 event at 0 ms.
        [EEG_erp, epochLabelReport] = set_LG_epoch_labels(EEG_erp, 'verbose', true); %#ok<NASGU>

        % NEW: deterministically exclude N+1 post-stimulation GS buffer epochs
        % from the primary EEG analysis. They remain present in the raw data.
        bufferMask = false(1,EEG_erp.trials);
        if isfield(EEG_erp.epoch,'eventlg_post_stim_buffer')
            bufferMask = logical([EEG_erp.epoch.eventlg_post_stim_buffer]);
        end
        nPostStimBufferExcluded = sum(bufferMask);
        if nPostStimBufferExcluded>0
            EEG_erp = pop_select(EEG_erp,'notrial',find(bufferMask));
            EEG_erp = eeg_checkset(EEG_erp);
            fprintf('   [%s %s] Excluded %d post-stimulation GS buffer epochs before EEG artifact rejection.\n', participant, thisExp, nPostStimBufferExcluded);
        end

        %% Apply the SAME session-level ICA rejections to this experiment
        if ~isempty(rowIdx) && iscell(T_LG.Bad_components_concatenated) ...
                && numel(T_LG.Bad_components_concatenated) >= rowIdx ...
                && ~isempty(T_LG.Bad_components_concatenated{rowIdx})

            EEG_erp = pop_subcomp(EEG_erp, T_LG.Bad_components_concatenated{rowIdx}, 0);
            EEG_erp = eeg_checkset(EEG_erp);
        end

        % (Optional) Inspect components on epoched data (comment in to review)
        % pop_viewprops(EEG_erp, 0, 1:size(EEG_erp.icawinv,2), 'freqrange', [2 80], 'ICLabel');

        %% Automatic epoch rejection options
        bad = [];

        % --- Peak-to-peak threshold (default in this pipeline)
        if ~isempty(autoRejPP) && autoRejPP > 0 && EEG_erp.trials>0
            p2p = zeros(EEG_erp.trials,1);
            for ep = 1:EEG_erp.trials
                seg = double(EEG_erp.data(:,:,ep));
                p2p(ep) = max(max(seg,[],2)-min(seg,[],2));
            end
            bad = find(p2p > autoRejPP);
        end

        % % --- Joint probability (comment-in to use; matches your earlier style)
        % threshJP = 3; % SD
        % [EEG_erp, locthresh, globthresh, idxJP] = pop_jointprob(EEG_erp, 1, 1:EEG_erp.nbchan, 3, 3, 0, 0, 0, [], 0);
        % bad = union(bad, find(EEG_erp.reject.rejjp));

        % % --- Kurtosis (comment-in to use)
        % [EEG_erp, ~, ~, idxK] = pop_rejkurt(EEG_erp, 1, 1:EEG_erp.nbchan, 3, 3, 0, 1, 0, [], 0);
        % bad = union(bad, find(EEG_erp.reject.rejkurt));

        % %--- Threshold (comment-in to use; example ±100 μV between 0–600 ms)
        % [EEG_erp, Indexes] = pop_eegthresh(EEG_erp, 1, 1:EEG_erp.nbchan, -40, 40, 0.4, 0.8, 1, 1);
        % bad = union(bad, find(EEG_erp.reject.rejthresh));

        % Apply epoch rejection if any flagged
        if ~isempty(bad)
            EEG_erp = pop_rejepoch(EEG_erp, bad, 0);
            EEG_erp = eeg_checkset(EEG_erp);
            fprintf('   [%s %s] Rejected %d epochs (auto rules)\n', participant, thisExp, numel(bad));
        end
        nKept = EEG_erp.trials;


        %% Interpolation (put back channels removed before ICA, if metadata exists)
        if isfield(EEG_erp.chaninfo,'removedchans') && ~isempty(EEG_erp.chaninfo.removedchans)
            EEG_erp = pop_interp(EEG_erp, EEG_erp.chaninfo.removedchans, 'spherical');
            EEG_erp = eeg_checkset(EEG_erp);
        end

        %% Re-reference (TP9/TP10 or AVG)
        if re_reference == true
            idTP9  = find(strcmpi({EEG_erp.chanlocs.labels}, 'TP9'), 1);
            idTP10 = find(strcmpi({EEG_erp.chanlocs.labels}, 'TP10'),1);
            EEG_erp = pop_reref(EEG_erp, [idTP9 idTP10]);    % Mastoids
            refLabel = 'TP9/TP10';
        else
            EEG_erp = pop_reref(EEG_erp, []);                % Average
            refLabel = 'AVG';
        end
        EEG_erp = eeg_checkset(EEG_erp);
        fprintf('   [%s %s] Re-referenced to %s\n', participant, thisExp, refLabel);

        %% Save ALL / GLOBAL STANDARD / GLOBAL DEVIANT datasets
        % Re-assert one scalar condition per epoch after all rejection steps.
        % This prevents neighboring T5 events inside the long epoch from being
        % mistaken for the trial-defining condition downstream.
        [EEG_erp, epochLabelReport] = set_LG_epoch_labels(EEG_erp, 'verbose', false); %#ok<NASGU>

        % Keep the historical filenames STANDARD/DEVIANT for compatibility,
        % but here they explicitly mean GLOBAL Standard/Deviant:
        %   GS = FLS + FLD
        %   GD = RLS + RLD
        taskCondition="";
        try
            [~,rawStem]=fileparts(dVhdr(1).name);
            taskMeta=read_eeg_visit_metadata(subFolder,rawStem,participant,analysisVisit,false);
            taskCondition=taskMeta.condition;
            EEG_erp.etc.erta_visit_metadata=taskMeta;
        catch ME
            warning('B3:TaskMetadata','New condition-specific figure unavailable: %s',ME.message);
        end

        taskAux=retain_aux_epochs(taskAux,EEG_erp);
        EEG_erp.etc.erta_aux=taskAux;
        tag = [participant '_' thisExp];
        EEG_erp.setname = ['auditory_all_LG_' tag];
        pop_saveset(EEG_erp, 'filename', ['auditory_all_LG_' tag '.set'], 'filepath', out_subfold);

        % Select by the scalar epoch label, NOT by "any event type in epoch".
        % The latter is unsafe for these long epochs because a neighboring
        % pattern can also fall inside the epoch window.
        epochEffect = string({EEG_erp.epoch.eventlg_effect});
        idxGS = find(epochEffect == "Global Standard");
        idxGD = find(epochEffect == "Global Deviant");

        EEGstd = pop_select(EEG_erp, 'trial', idxGS);
        EEGstd.etc.erta_aux=retain_aux_epochs(taskAux,EEGstd);
        EEGstd.setname = ['auditory_STANDARD_LG_' tag];
        pop_saveset(EEGstd, 'filename', ['auditory_STANDARD_LG_' tag '.set'], 'filepath', out_subfold);

        EEGdev = pop_select(EEG_erp, 'trial', idxGD);
        EEGdev.etc.erta_aux=retain_aux_epochs(taskAux,EEGdev);
        EEGdev.setname = ['auditory_DEVIANT_LG_' tag];
        pop_saveset(EEGdev, 'filename', ['auditory_DEVIANT_LG_' tag '.set'], 'filepath', out_subfold);

        nStd = EEGstd.trials;
        nDev = EEGdev.trials;

        %% Two channel pairs, each with old and new figure logic
        for pairNo=1:size(channelPairs,1)
        chans2plot=channelPairs(pairNo,:);
        pairSuffix='';
        if pairNo>1, pairSuffix=['_' upper(chans2plot{1}) '_' upper(chans2plot{2})]; end
        pairSuffix=regexprep(pairSuffix,'[^A-Za-z0-9_]','');
        % --- Bekinschtein (Local-global contrast)
        opts = struct();
        opts.toneOnsetsMs=toneOnsets_ms;
        opts.auxChannel=auxChannel;
        opts.auxChannelAliases=auxChannelAliases;
        opts.auxPreprocessing=auxPreprocessing;
        opts.auxRequired=auxRequired;
        opts.auxLineWidth=auxLineWidth;
        opts.auxControlColor=auxControlColor;
        opts.auxPairedColor=auxPairedColor;
        opts.expandForAux=expandPlotLimitsForAux;
        opts.topChan     = chans2plot{1};
        opts.botChan     = chans2plot{2};
        opts.winLocal    = topoLocalWin_ms;
        opts.winGlobal   = topoGlobalWin_ms;
        opts.tlim        = [epochWin*1000]; % x-axis range in ms
        opts.ylim        = erpPlotLimits;     % set [] to let MATLAB autoscale
        opts.topoGrid    = 80;         % 64–100 = crisp, older-EEGLAB look
        opts.topoShading = 'interp';   % no interpolation smoothing
        opts.topoContours= 8;
        opts.caxis = topoPlotLimits;       % colorbar scale

        opts.outPng      = fullfile(out_subfold, sprintf('Bekinschtein_%s_%s%s.png',participant,thisExp,pairSuffix));

        plot_bekinschtein_with_tones(EEG_erp, participant, thisExp, opts);
        close(gcf);
        drawnow;

        % New seven-curve ERP view and pure/paired/control topographies.
        if strlength(taskCondition)>0
            try
                newOpts=opts; newOpts.insetTopGap=topoInsetTopGap;
                newOpts.outPng=fullfile(out_subfold,sprintf('%s_V%d_%s_EEG_Pure_plus_%s_LocalGlobal%s.png', ...
                    participant,analysisVisit,thisExp,taskCondition,pairSuffix));
                [newFig,figureSummary]=plot_erta_eeg_control(EEG_erp,EEGcontrol,participant, ...
                    analysisVisit,thisExp,taskCondition,newOpts);
                save(strrep(newOpts.outPng,'.png','_summary.mat'),'figureSummary');
                close(newFig);
            catch ME
                warning('B3:NewFigure','%s %s new ERP figure failed: %s',participant,thisExp,ME.message);
            end
        end

        end % channel pairs: four requested figures per experiment

        % Keep the existing Mazancieux output on the first pair.
        opts.topChan=channelPairs{1,1}; opts.botChan=channelPairs{1,2};
        %%  Mazancieux ERPs contrast (local patters)
        % --- Mazancieux (two panels, four traces + SEM per panel)
        opts.ylim  = [-8 8];        % optional, to lock y-scale
        opts.outPng  = fullfile(out_subfold, sprintf('Mazancieux_%s_%s.png',participant,thisExp));
        plot_mz4(EEG_erp, participant, thisExp, opts);
        close(gcf);
        drawnow;

        %% Bookkeeping (per experiment)
        if ~isempty(rowIdx)
            if eidx == 1
                T_LG.Bad_trials_E1{rowIdx}       = bad;
                T_LG.Num_bad_trials_E1(rowIdx)   = numel(bad);
                T_LG.Num_good_trials_E1(rowIdx)  = nKept;
                T_LG.nStd_E1(rowIdx)             = nStd;
                T_LG.nDev_E1(rowIdx)             = nDev;
            else
                T_LG.Bad_trials_E2{rowIdx}       = bad;
                T_LG.Num_bad_trials_E2(rowIdx)   = numel(bad);
                T_LG.Num_good_trials_E2(rowIdx)  = nKept;
                T_LG.nStd_E2(rowIdx)             = nStd;
                T_LG.nDev_E2(rowIdx)             = nDev;
            end
        end

        % Save after each experiment (or move this save to after E1+E2)
        save(fullfile(outpath,'Info_EEG_LG.mat'), 'T_LG');
        fprintf('   [%s %s] kept: %d (Global Std=%d, Global Dev=%d)\n', participant, thisExp, nKept, nStd, nDev);

        % Write Excel export 
        excelFile = fullfile(outpath,'Info_EEG_LG.xlsx');
        write_TLG_excel_exact(T_LG, excelFile);
        fprintf('Excel export written: %s\n', excelFile);

        close all; drawnow;

    end % E1/E2 loop

    % Free memory before next subject
    clear EEG EEGc EEG_erp EEGstd EEGdev probs Plot_IC Brain Muscle Eye Heart LineN ChanN Other;
    drawnow;

    fprintf('   [%s] Pre-processing done!  \n', participant);

end % Participant loop

%%


function d = prefer_exact_stem(d, exactName)
if isempty(d), return; end
idx=find(strcmpi({d.name},exactName),1);
if ~isempty(idx), d=d([idx setdiff(1:numel(d),idx,'stable')]); end
end

function [EEGc,rej,lbltxt]=classify_eeg_components(EEGc,thr,out_subfold,figureTag)
    % ICLabel on concatenated data (session-level rejection list)
    lbltxt = strings(0);   % make safe even if ICLabel fails
    rej = [];

    try
        EEGc = iclabel(EEGc); % ICLabel
        % IC label components figure (Optional). Open IC property browser (comment in if you want to visually confirm):
        % pop_viewprops(EEGc, 0, 1:size(EEGc.icawinv,2), 'freqrange', [2 80], 'ICLabel'); % Probably one has to open it manually the first time to get the plugin installed
        % saveas(gcf , [out_subfold, 'All_ICA_', participant, '.jpg']); % Saving figure components

        % Pull probabilities & class names
        probs   = EEGc.etc.ic_classification.ICLabel.classifications;   % (IC x 7 classes)
        classes = EEGc.etc.ic_classification.ICLabel.classes;           % cellstr
        nIC     = size(EEGc.icawinv, 2);
        Plot_IC       = probs; % short alias for plotting

        % Any component with at least num% confidence of being an artifact or less than num% confidence if its brain
        Brain = probs(:,1);
        Muscle = probs(:,2);
        Eye = probs(:,3);
        Heart = probs(:,4);
        LineN = probs(:,5);
        ChanN = probs(:,6);
        Other = probs(:,7);

        rej = find( Eye>=thr.eye | Muscle>=thr.muscle | Heart>=thr.heart | LineN>=thr.line | ChanN>=thr.channel | Other>=thr.other | Brain<=thr.brainMin );
        rej = unique(rej(:)');

        % ===== Quick grid summary (scalp maps + dominant label) =====
        cols = 6; rows = ceil(nIC/cols);
        figGrid = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1], 'Renderer','opengl');
        tlg = tiledlayout(figGrid, rows, cols, 'Padding','compact', 'TileSpacing','compact');

        for k = 1:nIC
            nexttile;
            try
                topoplot(EEGc.icawinv(:,k), EEGc.chanlocs, 'numcontour', 4, 'electrodes', 'off');
            catch ME
                pop_topoplot(EEGc, 0, k, '', [1 1], 0, 'electrodes','off');
            end
            [pmax, imax] = max(Plot_IC(k,:));
            title(sprintf('IC%02d %s\n%.0f%%', k, classes{imax}, 100*pmax), 'FontSize', 8);
        end
        drawnow; pause(0.01);
        outpng = fullfile(out_subfold, ['All_ICA_grid_' figureTag '.png']);
        try
            exportgraphics(figGrid, outpng, 'Resolution', 200);
        catch ME
            saveas(figGrid, outpng);
        end

        close(figGrid);   % IMPORTANT: prevent figure/memory buildup
        drawnow;


        % Label list for bookkeeping
        classes = EEGc.etc.ic_classification.ICLabel.classes;
        lbltxt  = strings(size(rej));
        for i=1:numel(rej)
            [~, idxMax] = max(probs(rej(i),:));
            lbltxt(i) = string(classes{idxMax});
        end
        fprintf('   Session-level ICs flagged: %s\n', mat2str(rej));

        % (Optional) Save per-IC figures for flagged comps
        % bad_compons = rej;
        % for compy = 1:length(bad_compons)
        %     fig_IC = pop_prop_extended(EEGc, 0, bad_compons(compy)); %#ok<NASGU>
        %     saveas(gcf , fullfile(out_subfold, sprintf('Bad_compon_%d_%s_%s.png', ...
        %                   bad_compons(compy), char(lbltxt(compy)), participant)));
        %     close(gcf);
        % end

    catch ME
        warning('B3:ICLabelFailed','ICLabel failed for %s; no ICs rejected: %s',figureTag,ME.message);
        rej = [];
        lbltxt = strings(0);  % <-- add this line
    end


end

function EEGout=preprocess_control_eeg(dataRoot,outFolder,participant,visitNo,mainpath,HPF,LPF,epochWin,baseWin,thr,autoRejPP,re_reference,auxChannel,auxScale,auxAliases,auxMode,auxRequired)
EEGout=[];
folder=fullfile(dataRoot,participant,sprintf('Visit%d',visitNo));
d=dir(fullfile(folder,sprintf('%sV%dC*.vhdr',participant,visitNo)));
exact=strcmpi({d.name},sprintf('%sV%dC.vhdr',participant,visitNo));
if any(exact), d=d(exact); end
assert(~isempty(d),'No same-visit SxxV#C.vhdr found in %s.',folder);
assert(numel(d)==1,'Multiple control runs; keep one exact SxxV#C stem.');
[~,stem]=fileparts(d.name);
meta=read_eeg_visit_metadata(folder,stem,participant,visitNo,true);
icaFile=fullfile(outFolder,['ica_' participant '_C.set']);
assert(isfile(icaFile),'Control ICA is missing. Run updated B2 first: %s',icaFile);
EEGc=pop_loadset(icaFile); EEGc=eeg_checkset(EEGc);
assert(isfield(EEGc.etc,'control_preparation'),'Control preparation metadata missing; rerun updated B2.');
prep=EEGc.etc.control_preparation;
assert(strcmpi(prep.source_vhdr,d.name),'Control ICA and raw control filenames differ; rerun B2.');
[EEGc,rej,lbltxt]=classify_eeg_components(EEGc,thr,outFolder,[participant '_C']);

EEG=pop_loadbv(folder,d.name,[],[]); EEG=eeg_checkset(EEG);
EEGauxRaw=extract_aux_channel(EEG,auxChannel,auxScale,auxAliases,auxRequired,fullfile(folder,d.name));
EEG.ref='FCz';
EEG.etc.initial_reference=struct('label','FCz','kind','hardware','source','BrainVision online reference');
EEG=pop_chanedit(EEG,'lookup',[mainpath 'plugins\dipfit5.4\standard_BESA\standard-10-5-cap385.elp']);
rawLabels=string({EEG.chanlocs.labels});
assert(isequal(upper(rawLabels),upper(string(prep.original_channel_labels))), ...
    'Control montage differs from the ICA preparation montage; rerun B2.');
if ~isempty(prep.bad_channel_labels)
    [present,idx]=ismember(upper(string(prep.bad_channel_labels)),upper(rawLabels));
    assert(all(present),'Prepared bad channels not found in raw control.');
    EEG=pop_select(EEG,'rmchannel',idx); EEG=eeg_checkset(EEG);
end
assert(isequal(upper(string({EEG.chanlocs.labels})),upper(string({EEGc.chanlocs.labels}))), ...
    'Control ICA channel labels/order do not match raw control.');
EEG.icaact=[]; EEG.icawinv=EEGc.icawinv; EEG.icasphere=EEGc.icasphere;
EEG.icaweights=EEGc.icaweights; EEG.icachansind=EEGc.icachansind;

% Only in the stimulation-only recording is S10 itself the onset anchor.
% In E1/E2 the existing S9/T5 anchor remains unchanged.
codes=nan(1,numel(EEG.event));
for k=1:numel(EEG.event)
    typ=EEG.event(k).type;
    if isnumeric(typ) && isscalar(typ)
        codes(k)=double(typ);
    else
        token=regexp(strtrim(char(string(typ))),'^[Ss]\s*(\d+)$','tokens','once');
        if ~isempty(token), codes(k)=str2double(token{1}); end
    end
end
anchors=find(codes==10);
assert(~isempty(anchors),'No S10 stimulation-onset markers found in control BrainVision recording.');
assert(~any(codes==9),'S9 tone markers found in control; refusing to treat a task file as control.');
lat=double([EEG.event(anchors).latency]);
assert(all(isfinite(lat)) && all(diff(lat)>0),'Nonchronological/duplicate S10 control markers.');
if isfinite(meta.loggedControlBursts)
    assert(numel(anchors)==meta.loggedControlBursts, ...
        'Control has %d S10 markers but MAT logs %d bursts; inspect matching files.',numel(anchors),meta.loggedControlBursts);
end
for k=1:numel(anchors)
    EEG.event(anchors(k)).type='CONTROL_STIM';
    EEG.event(anchors(k)).control_marker_index=k;
end
EEG=eeg_checkset(EEG,'eventconsistency');
assert(numel(baseWin)==2 && all(isfinite(baseWin)) && baseWin(1)<baseWin(2) && ...
    baseWin(1)>=1000*epochWin(1) && baseWin(2)<=0,'Invalid control baseline.');
EEG=tag_aux_anchors(EEG);
controlAux=prepare_aux_epochs(EEGauxRaw,EEG,HPF,LPF,epochWin,baseWin,{'CONTROL_STIM'},auxMode);
clear EEGauxRaw;
EEG=pop_eegfiltnew(EEG,'locutoff',HPF);
EEG=pop_eegfiltnew(EEG,'hicutoff',LPF); EEG=eeg_checkset(EEG);
EEG=pop_epoch(EEG,{'CONTROL_STIM'},epochWin,'epochinfo','yes');
assert(EEG.trials>0,'No complete control epochs remain.');
EEG=pop_rmbase(EEG,baseWin); EEG=eeg_checkset(EEG);
EEG=label_control_epochs(EEG);
beforeRejection=[EEG.epoch.control_marker_index];
if ~isempty(rej), EEG=pop_subcomp(EEG,rej,0); EEG=eeg_checkset(EEG); end
bad=[];
if ~isempty(autoRejPP) && autoRejPP>0
    p2p=zeros(EEG.trials,1);
    for ep=1:EEG.trials
        seg=double(EEG.data(:,:,ep));
        p2p(ep)=max(max(seg,[],2)-min(seg,[],2));
    end
    bad=find(p2p>autoRejPP);
end
rejectedMarkerIndices=beforeRejection(bad);
assert(numel(bad)<EEG.trials,'All control epochs exceeded the EEG rejection threshold.');
if ~isempty(bad), EEG=pop_rejepoch(EEG,bad,0); EEG=eeg_checkset(EEG); end
if isfield(EEG.chaninfo,'removedchans') && ~isempty(EEG.chaninfo.removedchans)
    EEG=pop_interp(EEG,EEG.chaninfo.removedchans,'spherical'); EEG=eeg_checkset(EEG);
end
if re_reference
    idTP9=find(strcmpi({EEG.chanlocs.labels},'TP9'),1);
    idTP10=find(strcmpi({EEG.chanlocs.labels},'TP10'),1);
    assert(~isempty(idTP9) && ~isempty(idTP10),'Control mastoid reference channels are missing.');
    EEG=pop_reref(EEG,[idTP9 idTP10]); refLabel='TP9/TP10';
else
    EEG=pop_reref(EEG,[]); refLabel='AVG';
end
EEG=eeg_checkset(EEG); EEG=label_control_epochs(EEG);
EEG.etc.erta_visit_metadata=meta;
EEG.etc.erta_aux=retain_aux_epochs(controlAux,EEG);
controlQC=struct('source_vhdr',fullfile(folder,d.name),'n_markers',numel(anchors), ...
    'n_epoched',numel(beforeRejection),'n_rejected',numel(bad),'n_retained',EEG.trials, ...
    'unepoched_marker_indices',setdiff(1:numel(anchors),beforeRejection), ...
    'rejected_marker_indices',rejectedMarkerIndices, ...
    'retained_marker_indices',[EEG.epoch.control_marker_index], ...
    'bad_components',rej,'component_labels',strjoin(lbltxt,', '), ...
    'bad_channel_labels',{prep.bad_channel_labels},'baseline_ms',baseWin,'epoch_s',epochWin, ...
    'highpass_hz',HPF,'lowpass_hz',LPF,'peak_to_peak_uV',autoRejPP, ...
    'reference',refLabel,'anchor','S10 intended stimulation onset; no latency shift');
EEG.etc.erta_control=controlQC;
EEG.setname=['control_all_' participant];
EEG=pop_saveset(EEG,'filename',['control_all_' participant '.set'],'filepath',outFolder);
save(fullfile(outFolder,sprintf('%s_V%d_C_EEG_QC.mat',participant,visitNo)),'controlQC','meta');
fprintf('   [%s V%d CONTROL %s] %d/%d epochs retained; %d rejected, %d unepoched.\n', ...
    participant,visitNo,meta.condition,EEG.trials,numel(anchors),numel(bad),numel(anchors)-numel(beforeRejection));
EEGout=EEG;
end

function EEG=label_control_epochs(EEG)
for ep=1:EEG.trials
    eventIdx=numeric_list(EEG.epoch(ep).event);
    offsets=numeric_list(EEG.epoch(ep).eventlatency);
    assert(numel(eventIdx)==numel(offsets),'Control epoch event/latency mismatch.');
    atZero=abs(offsets)<=500/EEG.srate+1e-6;
    hit=[];
    for j=find(atZero)
        ev=eventIdx(j);
        if strcmp(string(EEG.event(ev).type),'CONTROL_STIM'), hit(end+1)=ev; end %#ok<AGROW>
    end
    assert(numel(hit)==1,'Control epoch lacks one unique zero-time S10 anchor.');
    EEG.epoch(ep).control_marker_index=EEG.event(hit).control_marker_index;
end
end

function x=numeric_list(v)
if iscell(v)
    x=nan(1,numel(v));
    for k=1:numel(v), x(k)=str2double(string(v{k})); end
else
    x=double(v(:)');
end
end

function meta=read_eeg_visit_metadata(folder,stem,participant,visitNo,isControl)
% Read acquired metadata only. Never call assignment-generating functions.
values=strings(0,1); sources=strings(0,1); nLogged=NaN;
matFile=fullfile(folder,[stem '.mat']);
if isfile(matFile)
    S=load(matFile,'experimentInfo','events_control','controlSequenceInfo');
    if isfield(S,'experimentInfo')
        I=S.experimentInfo;
        if isfield(I,'participantID'), assert(strcmpi(string(I.participantID),participant),'MAT participant mismatch.'); end
        if isfield(I,'visitNumber'), assert(I.visitNumber==visitNo,'MAT visit mismatch.'); end
        if isfield(I,'fileStem'), assert(strcmpi(string(I.fileStem),stem),'MAT file stem mismatch.'); end
        if isControl && isfield(I,'experiment'), assert(strcmpi(string(I.experiment),'VNS_CONTROL'),'MAT is not a control recording.'); end
        if isfield(I,'sessionCondition'), values=[values; string(I.sessionCondition)]; end
    end
    if isControl && isfield(S,'events_control') && ~isempty(S.events_control)
        E=S.events_control; nLogged=numel(E);
        if isfield(E,'condition'), values=[values; string({E.condition})']; end
        if isfield(E,'visit'), assert(all([E.visit]==visitNo),'Control burst visit mismatch.'); end
        if isfield(E,'marker_code'), assert(all([E.marker_code]==10),'Unexpected logged control marker code.'); end
    end
    if isControl && isfield(S,'controlSequenceInfo') && isfield(S.controlSequenceInfo,'condition')
        values=[values; string(S.controlSequenceInfo.condition)];
    end
    sources(end+1)=string(matFile);
end
% Same-run ASC can supply/check the recorded condition if needed.
ascFile=fullfile(folder,[stem '.asc']);
if isfile(ascFile)
    raw=fileread(ascFile);
    messages=regexp(raw,'(?m)^MSG[^\r\n]*','match');
    for k=1:numel(messages)
        c=regexp(messages{k},'\bCONDITION_(ACTIVE|SHAM)\b','tokens','once');
        if ~isempty(c), values(end+1,1)=string(c{1}); end %#ok<AGROW>
        v=regexp(messages{k},'\bVISIT_(\d+)\b','tokens','once');
        if ~isempty(v), assert(str2double(v{1})==visitNo,'ASC visit mismatch.'); end
    end
    sources(end+1)=string(ascFile);
end
values=upper(strtrim(values)); values=values(~ismissing(values) & values~="");
values=unique(values);
assert(numel(values)==1 && ismember(values,["ACTIVE","SHAM"]), ...
    'Missing/conflicting ACTIVE/SHAM metadata in same-run MAT/ASC for %s.',stem);
meta=struct('condition',values,'participant',string(participant),'visit',visitNo, ...
    'source',strjoin(sources,'; '),'fileStem',string(stem),'loggedControlBursts',nLogged);
end

function [fig,summary]=plot_erta_eeg_control(EEG,C,participant,visitNo,thisExp,condition,opts)
% All scalp/auxiliary averages use retained epochs. Aux control is S10-locked;
% aux paired pools every retained paired condition and is S9/T5-locked.
required={'eventlg_sound','eventlg_effect','eventlg_stim_paired','eventlg_post_stim_buffer'};
assert(EEG.trials>0 && all(isfield(EEG.epoch,required)),'Missing retained EEG trial labels.');
local=strings(1,EEG.trials); globalStatus=local; stim=false(1,EEG.trials); buffer=stim;
for ep=1:EEG.trials
    local(ep)=scalar_epoch_text(EEG.epoch(ep).eventlg_sound);
    globalStatus(ep)=scalar_epoch_text(EEG.epoch(ep).eventlg_effect);
    stim(ep)=scalar_epoch_boolean(EEG.epoch(ep).eventlg_stim_paired);
    buffer(ep)=scalar_epoch_boolean(EEG.epoch(ep).eventlg_post_stim_buffer);
end
assert(all(ismember(local,["Local Standard","Local Deviant"])) && ...
    all(ismember(globalStatus,["Global Standard","Global Deviant"])),'Unrecognized task epoch labels.');
pure=~buffer & ~stim; paired=~buffer & stim;
stimName="Sham"; if condition=="ACTIVE", stimName="taVNS"; end
controlAvailable=~isempty(C) && C.trials>0;
if controlAvailable
    assert(isfield(C.etc,'erta_visit_metadata') && C.etc.erta_visit_metadata.condition==condition, ...
        'Control/task condition mismatch.');
    assert(C.etc.erta_visit_metadata.visit==visitNo && ...
        strcmpi(C.etc.erta_visit_metadata.participant,participant),'Control participant/visit mismatch.');
end
[pairedAux,pairedTimes,pairedMeta]=aux_curve(EEG,paired);
controlAux=[]; controlTimes=[]; controlMeta=struct('available',false,'reason','Control unavailable');
if controlAvailable, [controlAux,controlTimes,controlMeta]=aux_curve(C,true(1,C.trials)); end
fprintf('   Figure AUX: control=%d trials; paired=%d trials.\n',size(controlAux,1),size(pairedAux,1));
if opts.auxRequired && isempty(pairedAux) && any(paired)
    error('B3:AuxEpochsMissing','Paired trials exist but AUX data are missing. Rerun this replacement B3 from its beginning; do not reuse old EEG variables.');
end
if opts.auxRequired && controlAvailable && isempty(controlAux)
    error('B3:AuxEpochsMissing','Control EEG exists but AUX data are missing. Rerun this replacement B3 from its beginning.');
end
times=double(EEG.times(:)'); channels={opts.topChan,opts.botChan}; kinds={'Local','Global'};
wins={opts.winLocal,opts.winGlobal};
masks={local=="Local Standard",local=="Local Deviant";globalStatus=="Global Standard",globalStatus=="Global Deviant"};
colors=[0 .447 .741;.85 .325 .098;0 .447 .741;.85 .325 .098;0 0 0;opts.auxControlColor;opts.auxPairedColor];
styles={'-','-','--','--','-','-','-'}; widths=[1.8 1.8 1.8 1.8 3.2 opts.auxLineWidth opts.auxLineWidth];
maps=cell(2,3); curveData=cell(2,7); curveTimes=cell(2,7); labels=cell(2,1); n=zeros(2,7);
selected=zeros(1,2); selectedC=zeros(1,2);
auxValues=nan(2,2); auxSEM=auxValues; auxWindowN=zeros(2,2);
for row=1:2
    selected(row)=find_channel(EEG,channels{row});
    if row==1, short=["LS","LD"]; else, short=["GS","GD"]; end
    labels{row}=["Pure "+short(1),"Pure "+short(2),stimName+" + "+short(1),stimName+" + "+short(2), ...
        "Pure "+stimName,"Artifact taVNS","Artifact taVNS paired"];
    mm=[pure & masks{row,1};pure & masks{row,2};paired & masks{row,1};paired & masks{row,2}];
    for k=1:4
        n(row,k)=nnz(mm(k,:)); curveTimes{row,k}=times;
        curveData{row,k}=reshape(double(EEG.data(selected(row),:,mm(k,:))),EEG.pnts,[])';
    end
    window=times>=wins{row}(1) & times<=wins{row}(2);
    assert(any(window),'No samples in requested topography window.');
    maps{row,1}=difference_map(EEG,window,mm(2,:),mm(1,:));
    maps{row,2}=difference_map(EEG,window,mm(4,:),mm(3,:));
    if controlAvailable
        selectedC(row)=find_channel(C,channels{row});
        curveData{row,5}=reshape(double(C.data(selectedC(row),:,:)),C.pnts,[])';
        curveTimes{row,5}=double(C.times(:)'); n(row,5)=C.trials;
        controlWindow=C.times>=wins{row}(1) & C.times<=wins{row}(2);
        assert(any(controlWindow),'No samples in control topography window.');
        maps{row,3}=mean(double(C.data(:,controlWindow,:)),[2 3]);
    end
    curveData{row,6}=controlAux; curveTimes{row,6}=controlTimes; n(row,6)=size(controlAux,1);
    curveData{row,7}=pairedAux; curveTimes{row,7}=pairedTimes; n(row,7)=size(pairedAux,1);
    [auxValues(row,1),auxSEM(row,1),auxWindowN(row,1)]=aux_window_summary(controlAux,controlTimes,wins{row});
    [auxValues(row,2),auxSEM(row,2),auxWindowN(row,2)]=aux_window_summary(pairedAux,pairedTimes,wins{row});
end
if isempty(opts.caxis)
    allMaps=[]; for k=1:numel(maps), allMaps=[allMaps;maps{k}(:)]; end %#ok<AGROW>
    maxAbs=max(abs(allMaps(isfinite(allMaps)))); if isempty(maxAbs) || maxAbs==0, maxAbs=1; end
    crng=[-maxAbs maxAbs];
else, crng=opts.caxis;
end
erpLimits=opts.ylim; allBounds=[]; auxBounds=[];
for row=1:2
    for k=1:7
        M=curveData{row,k}; if isempty(M), continue; end
        mu=mean(M,1,'omitnan'); nn=sum(isfinite(M),1);
        se=std(M,0,1,'omitnan')./sqrt(nn); se(nn<2)=0;
        b=[mu-se mu+se]; allBounds=[allBounds b]; %#ok<AGROW>
        if k>=6, auxBounds=[auxBounds b]; end %#ok<AGROW>
    end
end
allBounds=allBounds(isfinite(allBounds)); auxBounds=auxBounds(isfinite(auxBounds));
assert(~isempty(allBounds),'No finite waveforms to plot.');
if isempty(erpLimits)
    erpLimits=[min(allBounds) max(allBounds)]; pad=max(.1,.06*diff(erpLimits)); erpLimits=erpLimits+[-pad pad];
elseif opts.expandForAux && ~isempty(auxBounds)
    extent=max(abs(auxBounds))*1.06;
    erpLimits=[min(erpLimits(1),-extent) max(erpLimits(2),extent)];
end
% Wider panels; only the THREE actual scalp maps have cards and window strips.
fig=figure('Color','w','Units','pixels','Position',[20 40 2600 1120],'InvertHardcopy','off');
panelPositions=[.035 .105 .455 .80;.535 .105 .455 .80];
windowColors=[.35 .65 .95;.40 .72 .42;.30 .30 .30];
insetHeight=.265; insetY=1-opts.insetTopGap-insetHeight;
assert(isfinite(insetY) && insetY>0 && opts.insetTopGap>=0,'Invalid inset gap.');
stripY=.018+(0:2)*.044;
clippedMeans=zeros(2,7);
for row=1:2
    p=panelPositions(row,:); ax=axes('Parent',fig,'Position',p); hold(ax,'on');
    xlim(ax,opts.tlim); ylim(ax,erpLimits);
    set(ax,'FontSize',11,'Clipping','on','Color','w','Layer','bottom');
    draw_tone_lines(ax,opts.toneOnsetsMs);
    yline(ax,0,'-','Color',[.4 .4 .4],'HandleVisibility','off');
    available=~cellfun(@isempty,maps(row,:));
    for k=1:3
        if ~available(k), continue; end
        w=[max(wins{row}(1),opts.tlim(1)) min(wins{row}(2),opts.tlim(2))];
        if w(1)>=w(2), continue; end
        yy=erpLimits(1)+(stripY(k)+[0 .031])*diff(erpLimits);
        patch(ax,[w fliplr(w)],[yy(1) yy(1) yy(2) yy(2)],windowColors(k,:), ...
            'FaceAlpha',.20,'EdgeColor',windowColors(k,:),'LineWidth',.7,'HandleVisibility','off');
    end
    hh=gobjects(1,7);
    for k=1:7
        if n(row,k)>0
            hh(k)=draw_eeg_sem(ax,curveTimes{row,k},curveData{row,k},labels{row}(k),colors(k,:),styles{k},widths(k));
            mu=mean(curveData{row,k},1,'omitnan');
            clippedMeans(row,k)=nnz(mu<erpLimits(1) | mu>erpLimits(2));
        else
            hh(k)=plot(ax,NaN,NaN,'Color',colors(k,:),'LineStyle',styles{k},'LineWidth',widths(k), ...
                'DisplayName',char(labels{row}(k)+" (unavailable)"));
        end
    end
    % Draw every line above all SEM fills, with AUX lines uppermost.
    for k=1:7, uistack(hh(k),'top'); end
    grid(ax,'on'); box(ax,'off');
    xlabel(ax,'Time from T5 / control S10 marker (ms)'); ylabel(ax,'\muV');
    title(ax,sprintf('%s | Visit %d | %s | %s | %s | %s',participant,visitNo,thisExp, ...
        kinds{row},condition,upper(EEG.chanlocs(selected(row)).labels)),'Interpreter','none','FontSize',12);
    legend(ax,hh,'Location','southwest','Box','off','Interpreter','none','FontSize',10);
    set(ax,'Position',p);
    titles={sprintf('Pure tones (Dev-Std)\n%g-%g ms',wins{row}), ...
        sprintf('%s paired (Dev-Std)\n%g-%g ms',stimName,wins{row}), ...
        sprintf('Pure %s\n%g-%g ms',stimName,wins{row})};
    for k=1:3
        card=[.018+(k-1)*.227 insetY .207 insetHeight];
        if available(k)
            xx=opts.tlim(1)+(card(1)+card(3)/2)*diff(opts.tlim);
            yy=erpLimits(1)+card(2)*diff(erpLimits);
            stripCenter=erpLimits(1)+(stripY(k)+.0155)*diff(erpLimits);
            plot(ax,[xx mean(wins{row})],[yy stripCenter],'--','Color',windowColors(k,:)*.8, ...
                'LineWidth',.7,'HandleVisibility','off');
        end
        cx=opts.tlim(1)+[card(1) card(1)+card(3)]*diff(opts.tlim);
        cy=erpLimits(1)+[card(2) card(2)+card(4)]*diff(erpLimits);
        patch(ax,[cx fliplr(cx)],[cy(1) cy(1) cy(2) cy(2)],[1 1 1], ...
            'EdgeColor',windowColors(k,:),'LineWidth',.9,'HandleVisibility','off');
        text(ax,card(1)+card(3)/2,card(2)+card(4)-.01,titles{k},'Units','normalized', ...
            'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8.5, ...
            'FontWeight','bold','Color',windowColors(k,:)*.8,'Interpreter','none','Clipping','on');
        mapRect=[card(1)+.008 card(2)+.025 .152 .175];
        mapPosition=[p(1)+mapRect(1)*p(3) p(2)+mapRect(2)*p(4) mapRect(3)*p(3) mapRect(4)*p(4)];
        mapAx=axes('Parent',fig,'Position',mapPosition); %#ok<LAXES>
        if ~available(k)
            axis(mapAx,'off'); text(mapAx,.5,.5,'No retained data','HorizontalAlignment','center','FontSize',8);
        else
            if k==3, locations=C.chanlocs; channel=selectedC(row); else, locations=EEG.chanlocs; channel=selected(row); end
            plot_scalp_map(mapAx,maps{row,k}(:),locations,channel,crng,[{opts.auxChannel} opts.auxChannelAliases]);
            cb=colorbar(mapAx,'eastoutside'); cb.FontSize=7; cb.Ticks=unique([crng(1) 0 crng(2)]);
            cb.Position=[p(1)+(card(1)+.174)*p(3) p(2)+(card(2)+.043)*p(4) .005*p(3) .143*p(4)];
            set(mapAx,'Position',mapPosition);
            text(ax,card(1)+.177,card(2)+.196,'\muV','Units','normalized', ...
                'HorizontalAlignment','center','FontSize',8,'Interpreter','tex');
        end
    end
end
heading=sprintf('%s | Visit %d | %s | Pure tones + %s paired + control + auxiliary | Mean +/- SEM', ...
    participant,visitNo,thisExp,stimName);
annotation(fig,'textbox',[.02 .942 .96 .04],'String',heading,'EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',15,'FontWeight','bold','Interpreter','none');
auxNote=sprintf('AUX mode: %s. Recorded reference; same retained trials; no ICA subtraction. T1-T5 lines apply to task trials only.',opts.auxPreprocessing);
annotation(fig,'textbox',[.025 .013 .95 .045],'String',auxNote, ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',10,'Interpreter','none');
summary=struct('participant',participant,'visit',visitNo,'task',thisExp,'condition',condition, ...
    'channels',{channels},'local_labels',labels{1},'global_labels',labels{2},'retained_counts',n, ...
    'control_available',controlAvailable,'buffers_excluded_in_plot',nnz(buffer), ...
    'local_window_ms',opts.winLocal,'global_window_ms',opts.winGlobal, ...
    'difference','Deviant minus Standard','tone_onsets_ms',opts.toneOnsetsMs, ...
    'map_limits_uV',crng,'erp_limits_uV',erpLimits,'clipped_mean_samples',clippedMeans, ...
    'aux_control',controlMeta,'aux_paired',pairedMeta,'aux_window_mean_uV',auxValues, ...
    'aux_window_sem_uV',auxSEM,'aux_window_n',auxWindowN, ...
    'aux_window_rows',{{'Local','Global'}},'aux_window_columns',{{'Control','Paired'}}, ...
    'layout','Local left / Global right; three scalp maps per panel; AUX curves only');
if any(clippedMeans(:)), warning('B3:PlotClipping','Some mean waveforms exceed the requested voltage limits; see figure summary or set erpPlotLimits=[].'); end
try, exportgraphics(fig,opts.outPng,'Resolution',180); catch, saveas(fig,opts.outPng); end
end


function v=scalar_epoch_text(v)
if iscell(v), assert(numel(v)==1,'Ambiguous epoch label.'); v=v{1}; end
v=string(v); assert(isscalar(v) && ~ismissing(v),'Missing/scalar epoch label required.');
end

function v=scalar_epoch_boolean(v)
if iscell(v), assert(numel(v)==1,'Ambiguous epoch flag.'); v=v{1}; end
assert((isnumeric(v) || islogical(v)) && isscalar(v) && isfinite(v) && ...
    (v==0 || v==1),'Missing/invalid epoch flag.'); v=logical(v);
end

function idx=find_channel(EEG,label)
idx=find(strcmpi({EEG.chanlocs.labels},label));
assert(numel(idx)==1,'Requested channel %s is missing or duplicated.',label);
end

function m=difference_map(EEG,window,deviant,standard)
if ~any(deviant) || ~any(standard), m=[]; return; end
m=mean(double(EEG.data(:,window,deviant)),[2 3])-mean(double(EEG.data(:,window,standard)),[2 3]);
end

function h=draw_eeg_sem(ax,t,M,label,color,style,width)
mu=mean(M,1,'omitnan'); n=sum(isfinite(M),1);
sem=std(M,0,1,'omitnan')./sqrt(n); sem(n<2)=NaN;
fill(ax,[t fliplr(t)],[mu+sem fliplr(mu-sem)],color,'FaceAlpha',.15, ...
    'EdgeColor','none','HandleVisibility','off','Clipping','on');
h=plot(ax,t,mu,'Color',color,'LineStyle',style,'LineWidth',width,'Clipping','on', ...
    'DisplayName',sprintf('%s (n=%d)',label,size(M,1)));
end


function A=extract_aux_channel(EEG,label,scale,aliases,required,sourceFile)
% Capture the imported recorded channel BEFORE bad-channel removal/ICA/ref.
A=[]; labels={EEG.chanlocs.labels};
[idx,matchedLabel]=resolve_aux_channel(labels,label,aliases);
fprintf('   AUX source: %s\n',sourceFile);
if isempty(idx)
    message=sprintf('Requested AUX "%s" not found. Imported channels: %s. Check auxChannel/auxChannelAliases at the top.', ...
        label,strjoin(labels,', '));
    if required, error('B3:AuxMissing','%s',message); end
    warning('B3:AuxMissing','%s',message); return;
end
assert(isscalar(scale) && isfinite(scale) && scale>0,'Auxiliary unit conversion must be positive.');
A=pop_select(EEG,'channel',idx);
A.data=double(A.data)*scale;
assert(any(isfinite(A.data(:))),'Selected auxiliary channel has no finite samples.');
A.icaweights=[]; A.icasphere=[]; A.icawinv=[]; A.icaact=[]; A.icachansind=[];
A.etc.erta_aux_label=matchedLabel;
A.etc.erta_aux_scale=scale;
A.etc.erta_aux_source=sourceFile;
fprintf('   AUX selected: channel %d = "%s"; %d samples; %.3f Hz; unit scale %g.\n', ...
    idx,matchedLabel,EEG.pnts,EEG.srate,scale);
end

function [idx,matchedLabel]=resolve_aux_channel(labels,requested,aliases)
% Named aliases only. Never guess a channel number, location, or any AUX input.
if ischar(aliases), aliases={aliases}; end
wanted=[{char(requested)} cellstr(aliases(:))'];
canon=@(x) regexprep(lower(strtrim(char(x))),'[^a-z0-9]','');
keys=cellfun(canon,labels,'UniformOutput',false);
wanted=cellfun(canon,wanted,'UniformOutput',false);
idx=find(ismember(keys,unique(wanted)));
assert(numel(idx)<=1,'Multiple channels match the configured AUX names. Set auxChannel to the intended electrode and auxChannelAliases={} to disambiguate.');
matchedLabel=''; if ~isempty(idx), matchedLabel=char(labels{idx}); end
end

function EEG=tag_aux_anchors(EEG)
for k=1:numel(EEG.event), EEG.event(k).erta_aux_anchor_id=k; end
end

function A=prepare_aux_epochs(A,EEG,HPF,LPF,epochWin,baseWin,labels,mode)
if isempty(A), return; end
% Same continuous event sequence and epoch boundaries as the scalp data.
A.event=EEG.event; A.urevent=EEG.urevent;
label=A.etc.erta_aux_label; scale=A.etc.erta_aux_scale;
sourceFile=A.etc.erta_aux_source;
if strcmpi(mode,'erp')
    A=pop_eegfiltnew(A,'locutoff',HPF);
    A=pop_eegfiltnew(A,'hicutoff',LPF);
    processing='ERP bandpass and baseline only; no ICA subtraction/interpolation';
elseif strcmpi(mode,'raw')
    HPF=[]; LPF=[]; baseWin=[];
    processing='Recorded samples, epoching only; no filter/baseline/ICA/interpolation';
else, error('Unknown auxiliary preprocessing mode.');
end
A=pop_epoch(A,labels,epochWin,'epochinfo','yes');
if ~isempty(baseWin), A=pop_rmbase(A,baseWin); end
A=eeg_checkset(A);
ids=aux_epoch_ids(A,labels);
V=reshape(double(A.data),A.pnts,A.trials)';
A=struct('label',label,'times',double(A.times(:)'), 'data',V, ...
    'anchor_ids',ids,'anchor_labels',{labels},'highpass_hz',HPF,'lowpass_hz',LPF, ...
    'baseline_ms',baseWin,'unit_scale',scale,'units','microvolts', ...
    'reference','As recorded; no scalp re-reference', ...
    'processing',processing,'mode',mode,'source_vhdr',sourceFile, ...
    'selection','Pending scalp-trial rejection matching');
fprintf('   AUX "%s": %d epoched trials; mode=%s; baseline=%s ms.\n',label,size(V,1),mode,mat2str(baseWin));
end

function A=retain_aux_epochs(A,EEG)
if isempty(A), return; end
ids=aux_epoch_ids(EEG,A.anchor_labels);
assert(numel(unique(A.anchor_ids))==numel(A.anchor_ids),'Duplicate auxiliary epoch anchors.');
[ok,idx]=ismember(ids,A.anchor_ids);
assert(all(ok),'Cannot match retained scalp epochs to raw auxiliary epochs; refusing a misaligned overlay.');
assert(numel(A.times)==EEG.pnts && max(abs(A.times-double(EEG.times(:)')))<1e-6, ...
    'Auxiliary and scalp time grids differ.');
A.data=A.data(idx,:); A.anchor_ids=ids;
A.selection='Exactly the retained scalp epochs, matched by unique zero-time event IDs';
fprintf('   AUX "%s": matched %d retained scalp epochs.\n',A.label,numel(ids));
end

function ids=aux_epoch_ids(EEG,labels)
ids=nan(1,EEG.trials);
for ep=1:EEG.trials
    ev=numeric_list(EEG.epoch(ep).event); dt=numeric_list(EEG.epoch(ep).eventlatency);
    assert(numel(ev)==numel(dt),'Epoch event indices and latencies disagree.');
    hit=[];
    for k=find(abs(dt)<=500/EEG.srate+1e-6)
        if ismember(string(EEG.event(ev(k)).type),string(labels)), hit(end+1)=ev(k); end %#ok<AGROW>
    end
    assert(isscalar(hit),'Expected one unique zero-time auxiliary alignment anchor.');
    ids(ep)=EEG.event(hit).erta_aux_anchor_id;
end
end

function draw_tone_lines(ax,toneTimes)
% These lines apply to task tones; control-only trials contain no tones.
for k=1:numel(toneTimes)
    if toneTimes(k)<ax.XLim(1) || toneTimes(k)>ax.XLim(2), continue; end
    c=[.60 .60 .60]; style=':'; if k==5, c=[.3 .3 .3]; style='-'; end
    xline(ax,toneTimes(k),style,sprintf('T%d',k),'Color',c,'LineWidth',1, ...
        'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left', ...
        'FontSize',8,'HandleVisibility','off');
end
end

function [M,t,meta]=aux_curve(EEG,mask)
M=[]; t=[]; meta=struct('available',false,'reason','Auxiliary channel unavailable');
if isempty(EEG) || ~isfield(EEG,'etc') || ~isfield(EEG.etc,'erta_aux') || isempty(EEG.etc.erta_aux), return; end
A=EEG.etc.erta_aux;
assert(size(A.data,1)==EEG.trials && numel(mask)==EEG.trials,'Auxiliary retained-trial count mismatch.');
M=A.data(mask,:); t=A.times; meta=rmfield(A,'data');
meta.available=~isempty(M); meta.n_trials=size(M,1);
meta.anchor_ids=A.anchor_ids(mask);
if isempty(M), meta.reason='No retained trials for this auxiliary average'; end
end

function [value,se,n]=aux_window_summary(M,t,w)
value=NaN; se=NaN; n=0;
if isempty(M), return; end
mask=t>=w(1) & t<=w(2); if ~any(mask), return; end
v=mean(M(:,mask),2,'omitnan'); v=v(isfinite(v)); n=numel(v);
if n>0, value=mean(v); end
if n>1, se=std(v,0)/sqrt(n); end
end

function plot_scalp_map(ax,values,locations,selected,limits,auxLabel)
% The single ear auxiliary channel is not a scalp sample, even if a location
% was accidentally assigned by a lookup. Never use it to interpolate a map.
keys=regexprep(lower(string({locations.labels})),'[^a-z0-9]','');
excluded=regexprep(lower(string(auxLabel)),'[^a-z0-9]','');
keep=~ismember(keys,excluded);
assert(keep(selected),'Selected ERP channel must be a scalp channel.');
newSelected=sum(keep(1:selected));
axes(ax);
topoplot(double(values(keep)),locations(keep),'electrodes','on', ...
    'emarker2',{newSelected,'o',[0 0 0],5,1},'numcontour',6,'gridscale',120, ...
    'maplimits',limits,'whitebk','on');
caxis(ax,limits);
end


function plot_bekinschtein_with_tones(EEG, participant, thisExp, opts)
% Local (MMN) @ topChan vs Global (P3b) @ botChan
% Uses EEG.epoch.eventlg_sound / eventlg_effect (== after your rejections)

%% ---------- defaults (NO anonymous helpers) ----------
if nargin < 4, opts = struct; end

if ~isfield(opts,'topChan')   || isempty(opts.topChan),    opts.topChan   = 'Fz';   end
if ~isfield(opts,'botChan')   || isempty(opts.botChan),    opts.botChan   = 'Pz';   end
if ~isfield(opts,'winLocal')  || isempty(opts.winLocal),   opts.winLocal  = [130 200]; end
if ~isfield(opts,'winGlobal') || isempty(opts.winGlobal),  opts.winGlobal = [320 600]; end
if ~isfield(opts,'outPng'),                                opts.outPng    = '';        end
if ~isfield(opts,'tlim')      || isempty(opts.tlim),       opts.tlim      = [EEG.times(1) EEG.times(end)]; end
if ~isfield(opts,'ylim'),                                 opts.ylim      = [];        end

times = EEG.times;

%% ---------- masks ----------
[mLocStd, mLocDev] = masks_from_epoch_field_local(EEG);
[mGloStd, mGloDev] = masks_from_epoch_field_global(EEG);

nLocStd = nnz(mLocStd); nLocDev = nnz(mLocDev);
nGloStd = nnz(mGloStd); nGloDev = nnz(mGloDev);

%% ---------- channels ----------
% logical masks
ciTop = strcmpi({EEG.chanlocs.labels}, opts.topChan);
ciBot = strcmpi({EEG.chanlocs.labels}, opts.botChan);

assert(nnz(ciTop)==1 && nnz(ciBot)==1,'Requested scalp channel missing/duplicated; no fallback channel is used.');

% convert logical -> numeric index (for markers)
idxTop = find(ciTop,1);
idxBot = find(ciBot,1);

%% ---------- ERPs ----------
Lstd = mean(double(EEG.data(idxTop,:,mLocStd)), 3);
Ldev = mean(double(EEG.data(idxTop,:,mLocDev)), 3);
Gstd = mean(double(EEG.data(idxBot,:,mGloStd)), 3);
Gdev = mean(double(EEG.data(idxBot,:,mGloDev)), 3);

%% ---------- topography data ----------
iL = times >= opts.winLocal(1)  & times <= opts.winLocal(2);
iG = times >= opts.winGlobal(1) & times <= opts.winGlobal(2);

diffLocal  = mean(double(EEG.data(:,iL,mLocDev)), [2 3]) - ...
    mean(double(EEG.data(:,iL,mLocStd)), [2 3]);
diffGlobal = mean(double(EEG.data(:,iG,mGloDev)), [2 3]) - ...
    mean(double(EEG.data(:,iG,mGloStd)), [2 3]);


% After diffLocal / diffGlobal are computed --> Color scale for topoplot
if isfield(opts,'caxis') && ~isempty(opts.caxis)
    crng = opts.caxis;                        % use caller-provided range
else
    cmax = max(abs([diffLocal; diffGlobal])); % default: from the two windows
    if cmax==0, cmax=1; end
    crng = [-cmax cmax];
end


%% ---------- figure ----------
fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w');
t   = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

%% ====== 1) Local trace ======
ax1 = nexttile(t,1); hold(ax1,'on');
ax1.FontSize = 11;   % tick numbers on both axes for the top ERP
set(groot,'defaultAxesFontSize',11);  % affects all future axes

p1 = plot(times, Lstd, 'LineWidth', 2, 'Color',[0 0.45 0.85]); %blue
p2 = plot(times, Ldev, 'LineWidth', 2, 'Color', [0.90 0.45 0.10]); %orange
yline(ax1,0,'k-');

if ~isempty(opts.ylim), ylim(ax1, opts.ylim); end
xlim(ax1, opts.tlim);

yl = ylim(ax1);
% shaded MMN window (clipped to tlim)
winL_clip = [max(opts.winLocal(1), opts.tlim(1)) , ...
    min(opts.winLocal(2), opts.tlim(2))];
if winL_clip(1) < winL_clip(2)
    patch(ax1, [winL_clip fliplr(winL_clip)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.82 0.90 1.00], 'FaceAlpha',0.5, 'EdgeColor','none','HandleVisibility','off');
end

draw_tone_lines(ax1,opts.toneOnsetsMs);
grid(ax1,'on');
xlabel(ax1,'Time (ms)', FontSize=12);
ylabel(ax1,'\muV', FontSize=12);

% Local title
if isfield(EEG,'etc') && isfield(EEG.etc,'grandavg') && isfield(EEG.etc.grandavg,'nSubTotal')
    title(ax1, sprintf('%s | Local (MMN expected)   [Nsub=%d | Ntrials Std=%d Dev=%d]', ...
        upper(EEG.chanlocs(idxTop).labels), ...
        EEG.etc.grandavg.nSubTotal, ...
        EEG.etc.grandavg.nTrialLocalStd, EEG.etc.grandavg.nTrialLocalDev), ...
        'Interpreter','none', FontSize=12);
else
    title(ax1, sprintf('%s | Local (MMN expected)   [Std=%d, Dev=%d]', ...
        upper(EEG.chanlocs(idxTop).labels), nLocStd, nLocDev), ...
        'Interpreter','none', FontSize=12);
end


legend(ax1, [p1 p2], {'Local Standard','Local Deviant'}, ...
    'Location','southwest','Box','off', FontSize=12);

%% ====== 2) Local topo ======
ax2 = nexttile(t,2); cla(ax2); axes(ax2);
idxTop = find(ciTop,1);                        % numeric index

% optional: silence position warning from topoplot in tiledlayout
ws = warning('off','MATLAB:ui:Axes:UnsupportedConveniencePositionSet');

plot_scalp_map(ax2,diffLocal(:),EEG.chanlocs,idxTop,crng,[{opts.auxChannel} opts.auxChannelAliases]);

warning(ws);  % restore

% derive the shown window from the logical mask (no reliance on winLocal var)
wL = [min(times(iL)) max(times(iL))];
title(sprintf('Local at %d–%d ms', round(wL(1)), round(wL(2))), FontSize=12);
cb = colorbar(ax2,'eastoutside'); cb.Label.String = '\muV (Dev-Std)', cb.Label.FontSize = 12;

% colorbar scale topoplots:
caxis(ax2, crng);


%% ====== 3) Global trace ======
ax3 = nexttile(t,3); hold(ax3,'on');
ax3.FontSize = 11;   % bottom ERP
set(groot,'defaultAxesFontSize',11);  % affects all future axes

p3 = plot(times, Gstd, 'LineWidth', 2, 'Color',[0 0.45 0.85]);
p4 = plot(times, Gdev, 'LineWidth', 2, 'Color',[0.90 0.45 0.10]);
yline(ax3,0,'k-');

if ~isempty(opts.ylim), ylim(ax3, opts.ylim); end
xlim(ax3, opts.tlim);

yl = ylim(ax3);
winG_clip = [max(opts.winGlobal(1), opts.tlim(1)) , ...
    min(opts.winGlobal(2), opts.tlim(2))];
if winG_clip(1) < winG_clip(2)
    patch(ax3, [winG_clip fliplr(winG_clip)], [yl(1) yl(1) yl(2) yl(2)], ...
        [1.00 0.90 0.82], 'FaceAlpha',0.5, 'EdgeColor','none','HandleVisibility','off');
end

draw_tone_lines(ax3,opts.toneOnsetsMs);
grid(ax3,'on');
xlabel(ax3,'Time (ms)', FontSize=12);
ylabel(ax3,'\muV', FontSize=12);



% Global title
if isfield(EEG,'etc') && isfield(EEG.etc,'grandavg') && isfield(EEG.etc.grandavg,'nSubTotal')
    title(ax3, sprintf('%s | Global (P3b expected)   [Nsub=%d | Ntrials Std=%d Dev=%d]', ...
        upper(EEG.chanlocs(idxBot).labels), ...
        EEG.etc.grandavg.nSubTotal, ...
        EEG.etc.grandavg.nTrialGlobalStd, EEG.etc.grandavg.nTrialGlobalDev), ...
        'Interpreter','none', FontSize=12);
else
    title(ax3, sprintf('%s | Global (P3b expected)   [Std=%d, Dev=%d]', ...
        upper(EEG.chanlocs(idxBot).labels), nGloStd, nGloDev), ...
        'Interpreter','none', FontSize=12);
end

legend(ax3, [p3 p4], {'Global Standard','Global Deviant'}, ...
    'Location','southwest','Box','off', FontSize=12);


%% ====== 4) Global topo ======
ax4 = nexttile(t,4); cla(ax4); axes(ax4);
idxBot = find(ciBot,1);

ws = warning('off','MATLAB:ui:Axes:UnsupportedConveniencePositionSet');

plot_scalp_map(ax4,diffGlobal(:),EEG.chanlocs,idxBot,crng,[{opts.auxChannel} opts.auxChannelAliases]);

warning(ws);

wG = [min(times(iG)) max(times(iG))];
title(sprintf('Global at %d–%d ms', round(wG(1)), round(wG(2))), FontSize=12);
cb = colorbar(ax4,'eastoutside'); cb.Label.String = '\muV (Dev-Std)', cb.Label.FontSize = 12;
% colorbar scale topoplots:
caxis(ax4, crng);


%% ====== super title ======
sgtitle(sprintf('%s | %s', participant, thisExp), ...
    'Color', "k", 'FontWeight', 'bold', 'Interpreter', 'none', FontSize=14);

%% ====== turn off toolbars (kill warning) ======
axs = findall(fig,'Type','axes');
for k = 1:numel(axs)
    try
        axtoolbar(axs(k), 'none');
    catch
    end
end


%% ====== export ======
if ~isempty(opts.outPng)
    try
        exportgraphics(fig, opts.outPng, 'Resolution', 200);
    catch
        saveas(fig, opts.outPng);
    end
end

end

%% ===== helper: local masks (safe for Local only) =====
function [mStd, mDev] = masks_from_epoch_field_local(EEG)
mStd = false(1, EEG.trials);
mDev = false(1, EEG.trials);

% 1) Preferred: explicit local labels
if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_sound')
    vals = {EEG.epoch.eventlg_sound};
    for i = 1:numel(vals)
        v = vals{i};
        % normalize to a row of strings
        if ischar(v) || isstring(v)
            vv = string(v);
        elseif iscell(v)
            vv = string(v(:)');
        else
            vv = string.empty(1,0);
        end
        mStd(i) = any(strcmpi(vv,'Local Standard'));
        mDev(i) = any(strcmpi(vv,'Local Deviant'));
    end
    return;
end

% 2) Fallback: parse Mazancieux-style combined labels in eventlg_type
if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type')
    canon = @(s) lower(strtrim(regexprep(strrep(strrep(string(s),'pattern',''),'_',' '),'\s+',' ')));
    for i = 1:numel(EEG.epoch)
        ci = "";
        try
            ci = canon(EEG.epoch(i).eventlg_type);
        catch
            % leave as "", will evaluate to false below
        end
        mStd(i) = contains(ci, "local standard");
        mDev(i) = contains(ci, "local deviant");
    end
    return;
end

% 3) If we get here, we cannot infer local labels safely
error(['masks_from_epoch_field_local: No local labels found. ' ...
       'Need EEG.epoch.eventlg_sound or eventlg_type containing "local standard/deviant".']);
end

%% ===== helper: global masks =====
function [mStd, mDev] = masks_from_epoch_field_global(EEG)
mStd = false(1, EEG.trials);
mDev = false(1, EEG.trials);

if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_effect')
    vals = {EEG.epoch.eventlg_effect};
    for i = 1:numel(vals)
        v = vals{i};
        if ischar(v) || isstring(v)
            vv = string(v);
        elseif iscell(v)
            vv = string(v(:)');
        else
            vv = string.empty(1,0);
        end
        mStd(i) = any(strcmpi(vv, 'Global Standard'));
        mDev(i) = any(strcmpi(vv, 'Global Deviant'));
    end
else
    % fallback from eventtype
    et = {EEG.epoch.eventtype};
    for i = 1:numel(et)
        v = et{i};
        if ischar(v) || isstring(v)
            vv = string(v);
        elseif iscell(v)
            vv = string(v(:)');
        else
            vv = string.empty(1,0);
        end
        mStd(i) = any(strcmpi(vv,'STANDARD'));
        mDev(i) = any(strcmpi(vv,'DEVIANT'));
    end
end
end

%% ===== helper: STANDARD/DEVIANT fallback =====
function [mStd, mDev] = masks_from_eventtype_stddev(EEG)
mStd = false(1, EEG.trials);
mDev = false(1, EEG.trials);
if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventtype')
    vals = {EEG.epoch.eventtype};
    for i = 1:numel(vals)
        v = vals{i};
        if ischar(v) || isstring(v)
            mStd(i) = strcmpi(v,'STANDARD');
            mDev(i) = strcmpi(v,'DEVIANT');
        elseif iscell(v)
            mStd(i) = any(strcmpi(v,'STANDARD'));
            mDev(i) = any(strcmpi(v,'DEVIANT'));
        end
    end
end
end

%% ===== helper: pick-like =====
function mask = pick_like_local(allLabs, prefs)
mask = false(1,numel(allLabs));
AL = upper(string(allLabs));
for k = 1:numel(prefs)
    hit = find(AL == upper(string(prefs{k})), 1);
    if ~isempty(hit)
        mask(hit) = true;
        return;
    end
end
end
