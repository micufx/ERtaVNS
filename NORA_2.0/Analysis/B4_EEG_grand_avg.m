clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% Grand average (Local–Global | condition-aware, compatible with your plots)
% Uses auditory_all_LG_* so we keep 4-condition LG labels:
%   1) Frequent Local Standard
%   2) Rare     Local Standard
%   3) Frequent Local Deviant
%   4) Rare     Local Deviant
%
% Builds a grand dataset where each epoch = SUBJECT-LEVEL ERP for one of those 4 labels,
% and sets EEG.epoch.eventlg_type / eventlg_sound / eventlg_effect / eventtype.
%
% Also stores totals for plotting:
%   EEGga.etc.grandavg.nSubTotal
%   EEGga.etc.grandavg.nTrialLocalStd / nTrialLocalDev
%   EEGga.etc.grandavg.nTrialGlobalStd / nTrialGlobalDev
%
% FIX (Miguel): Stores correct subject_id per epoch by tracking subjects PER condition.
%               This avoids wrong mapping when some subjects have 0 trials in a condition.
%
% Miguel Contreras-Altamirano, 2026

%% Paths & files
mainpath = 'L:\Downloads\eeglab2023.0\';                      % EEGLAB folder
outpath  = 'L:\Downloads\ERtaVNS\NORA_2.0_PARAMETERIZED_READY\NORA_2.0\Analysis_Output';  % output root

skipSubjects_all = {}; % {'S02','S13'};

% Add EEGLAB properly
if exist('eeglab','file') ~= 2
    addpath(mainpath);
end

grandFold = fullfile(outpath, 'GRAND_AVG');
if ~exist(grandFold, 'dir'), mkdir(grandFold); end

% Experiments
expList = {'E1','E2'};

% Plot options
chans2plot = {'CZ','PZ'};

% The 4 LG condition labels expected by plot_mz4
targets = ["Frequent Local Standard","Rare Local Standard", ...
    "Frequent Local Deviant","Rare Local Deviant"];


%% ERP evolution
% --- Timtopo export options (NEW) ---
optsTT = struct();
optsTT.limits    = [-2500 700];                  % what timtopo shows on x-axis
optsTT.plottimes = linspace(-1800,-1300,13);     % default 13 topos (EDIT here)
optsTT.exportDPI = 300;
optsTT.doSubjects = true;  % per participant
optsTT.doGrand    = true;  % grand average

%% ===================== OPTIONAL SKIP LIST =====================
% Put subjects here if you want to exclude them from GRAND AVERAGE.
% Examples:
%   skipSubjects_all = {'S02','S13'};
%   skipSubjects_byExp.E1 = {'S02'};
%   skipSubjects_byExp.E2 = {'S07','S10'};
%
% Notes:
% - Matching is case-insensitive.
% - Skipping is applied BEFORE loading a subject.
% - Skipping cannot "shift rows" because we store included subjects by ID.

skipSubjects_byExp = struct();   % <-- FIX: don't create bogus fields
skipSubjects_byExp.E1 = {};      % only for E1
skipSubjects_byExp.E2 = {};      % only for E2

%% --- Find participants that have ALL sets (at least E1) ---
dSubs = dir(outpath);
isSub = [dSubs.isdir] & ~ismember({dSubs.name},{'.','..','GRAND_AVG'});
subNames = {dSubs(isSub).name};
subNames = subNames(~cellfun('isempty', regexp(subNames, '^S\d+', 'once'))); % keep S##

participant_list = {};
for si = 1:length(subNames)
    p = subNames{si};
    fAllE1 = fullfile(outpath, p, sprintf('auditory_all_LG_%s_E1.set', p));
    if exist(fAllE1,'file')
        participant_list{end+1} = p; %#ok<SAGROW>
    end
end

if isempty(participant_list)
    error('No participants found with expected: out\\Sxx\\auditory_all_LG_Sxx_E1.set');
end
participant_list = sort(participant_list);
fprintf('Participants found for grand average (>=E1 ALL): %d\n', numel(participant_list));

%% EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU> % Open eeglab

%% Loop experiments
for eidx = 1:numel(expList)

    thisExp = expList{eidx};
    fprintf('\n=== Grand average for %s ===\n', thisExp);

    %% Build active skip list (per exp)
    skipNow = string(skipSubjects_all);
    if isfield(skipSubjects_byExp, thisExp)
        skipNow = [skipNow, string(skipSubjects_byExp.(thisExp))]; %#ok<AGROW>
    end
    skipNow = unique(lower(strtrim(skipNow)));
    if ~isempty(skipNow)
        fprintf('Skip list for %s: %s\n', thisExp, strjoin(cellstr(skipNow), ', '));
    end

    %% Template
    templateEEG    = [];
    templateLabels = {};

    % Store subject-level ERPs per condition (each epoch = 1 subject ERP)
    condERP  = cell(1,4);              % [chan x pnts x nSubHavingThisCond]
    condSubj = cell(1,4);              % [nSubHavingThisCond x 1] subject IDs aligned to condERP
    for k = 1:4
        condERP{k}  = [];
        condSubj{k} = strings(0,1);
    end

    % Reset bookkeeping PER EXPERIMENT (trial counts are per included subject)
    trialCounts = zeros(0,4);          % rows=subjects included, cols=4 LG bins (FLS,RLS,FLD,RLD)
    subjIDs     = strings(0,1);

    %% Subject loop
    for si = 1:length(participant_list)

        participant = participant_list{si};

        % Skip check
        if ismember(lower(string(participant)), skipNow)
            fprintf('   [%s %s] In skip list -> skipping subject.\n', participant, thisExp);
            continue;
        end

        subFold = fullfile(outpath, participant);
        fAll    = fullfile(subFold, sprintf('auditory_all_LG_%s_%s.set', participant, thisExp));
        if ~exist(fAll,'file')
            fprintf('   [%s %s] Missing ALL set. Skipping.\n', participant, thisExp);
            continue;
        end


        EEGall = [];
        while isempty(EEGall)
            try
                EEGall = pop_loadset(fAll);
                EEGall = eeg_checkset(EEGall);
            catch
                warning('pop_loadset failed for %s. Retrying once...', fAll);
                pause(2);
                EEGall = pop_loadset(fAll);
                EEGall = eeg_checkset(EEGall);
            end
        end


        % Init template from first valid subject
        if isempty(templateEEG)
            templateEEG    = EEGall;
            templateLabels = {templateEEG.chanlocs.labels};
            fprintf('   Template set: %s (%d chans, %d pnts)\n', participant, templateEEG.nbchan, templateEEG.pnts);
        end

        % Reorder channels to template
        try
            EEGall = reorder_to_template(EEGall, templateLabels);
        catch
            warning('   [%s %s] Channel mismatch. Skipping.', participant, thisExp);
            continue;
        end

        % Safety checks
        if EEGall.pnts ~= templateEEG.pnts
            warning('   [%s %s] pnts mismatch. Skipping.', participant, thisExp);
            continue;
        end
        if numel(EEGall.times) ~= numel(templateEEG.times) || any(abs(EEGall.times - templateEEG.times) > 1e-6)
            warning('   [%s %s] times mismatch. Skipping.', participant, thisExp);
            continue;
        end

        % Need eventlg_type to build 4 condition masks
        if ~(isfield(EEGall,'epoch') && isfield(EEGall.epoch,'eventlg_type'))
            warning('   [%s %s] EEG.epoch.eventlg_type missing. Skipping.', participant, thisExp);
            continue;
        end

        % Build masks for 4 conditions
        M = false(4, EEGall.trials);
        for ep = 1:EEGall.trials
            lablist = to_strlist(EEGall.epoch(ep).eventlg_type);
            for k = 1:4
                M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
            end
        end

        % Count real trials per bin for THIS subject
        tFLS = nnz(M(1,:));
        tRLS = nnz(M(2,:));
        tFLD = nnz(M(3,:));
        tRLD = nnz(M(4,:));

        % Only include subject if they contribute something
        if (tFLS + tRLS + tFLD + tRLD) == 0
            warning('   [%s %s] No trials matched the 4 LG bins. Skipping subject.', participant, thisExp);
            continue;
        end

        % Save per-subject trial counts (for totals)
        trialCounts(end+1,:) = [tFLS tRLS tFLD tRLD]; %#ok<SAGROW>
        subjIDs(end+1,1)     = string(participant);   %#ok<SAGROW>

        % Compute subject ERP per condition; append ERP + subject ID PER condition
        for k = 1:4
            msk = M(k,:);
            if nnz(msk) == 0
                continue;
            end

            thisERP = mean(double(EEGall.data(:,:,msk)), 3); % [chan x pnts]

            if isempty(condERP{k})
                condERP{k} = zeros(size(thisERP,1), size(thisERP,2), 0);
            end

            condERP{k}(:,:,end+1)  = thisERP;              %#ok<SAGROW>
            condSubj{k}(end+1,1)   = string(participant);  %#ok<SAGROW>
        end

        fprintf('   [%s %s] Included trials: FLS=%d RLS=%d FLD=%d RLD=%d\n', ...
            participant, thisExp, tFLS, tRLS, tFLD, tRLD);

    end % subject loop

    %% Trial totals for plotting in Bekinschtein
    if isempty(subjIDs)
        warning('No valid subjects for %s. Skipping.', thisExp);
        continue;
    end

    nSubTotal = numel(subjIDs);

    % Local: Std=(FLS+RLS), Dev=(FLD+RLD)
    nTrialLocalStd = sum(trialCounts(:,1) + trialCounts(:,2));
    nTrialLocalDev = sum(trialCounts(:,3) + trialCounts(:,4));

    % Global (your mapping): Frequent = Global Standard, Rare = Global Deviant
    % Std=(FLS+FLD), Dev=(RLS+RLD)
    nTrialGlobalStd = sum(trialCounts(:,1) + trialCounts(:,3));
    nTrialGlobalDev = sum(trialCounts(:,2) + trialCounts(:,4));

    %% Build GRAND dataset (epochs = subject ERPs per LG condition)
    nbchan = templateEEG.nbchan;
    pnts   = templateEEG.pnts;

    nPerCond = zeros(1,4);
    for k = 1:4
        if ~isempty(condERP{k}), nPerCond(k) = size(condERP{k},3); end
    end
    nTotal = sum(nPerCond);

    if nTotal == 0
        warning('No valid subject ERPs for %s. Skipping.', thisExp);
        continue;
    end

    fprintf('-> %s epochs: total=%d | FLS=%d RLS=%d FLD=%d RLD=%d\n', ...
        thisExp, nTotal, nPerCond(1), nPerCond(2), nPerCond(3), nPerCond(4));

    dataGA = zeros(nbchan, pnts, nTotal);

    % Build epoch structs WITH the required fields (+ subject_id)
    epoch = repmat(struct('eventlg_type','','eventlg_sound','','eventlg_effect','', ...
        'eventtype','', 'subject_id',''), 1, nTotal);

    % Map each k to local/global strings for Bekinschtein masks
    localSound = ["Local Standard","Local Standard","Local Deviant","Local Deviant"];
    globalEff  = ["Global Standard","Global Deviant","Global Standard","Global Deviant"];

    epIdx = 0;
    for k = 1:4
        if isempty(condERP{k}), continue; end
        nK = size(condERP{k},3);

        % Sanity: subject list must match ERP stack depth
        if numel(condSubj{k}) ~= nK
            warning('[%s %s] condSubj{%d} length (%d) != nK (%d). Fixing by truncation.', ...
                'GRAND_AVG', thisExp, k, numel(condSubj{k}), nK);
            condSubj{k} = condSubj{k}(1:min(end,nK));
            nK = min(nK, numel(condSubj{k}));
        end

        for i = 1:nK
            epIdx = epIdx + 1;
            dataGA(:,:,epIdx) = condERP{k}(:,:,i);

            epoch(epIdx).eventlg_type   = char(targets(k));
            epoch(epIdx).eventlg_sound  = char(localSound(k));
            epoch(epIdx).eventlg_effect = char(globalEff(k));

            % Correct subject ID for THIS condition-epoch
            epoch(epIdx).subject_id     = char(condSubj{k}(i));

            % eventtype fallback (for STANDARD/DEVIANT-only code paths)
            if k==1 || k==2
                epoch(epIdx).eventtype = 'STANDARD';
            else
                epoch(epIdx).eventtype = 'DEVIANT';
            end
        end
    end

    % Final sanity
    if epIdx ~= nTotal
        warning('[%s %s] Built epochs epIdx=%d but expected nTotal=%d. Truncating to epIdx.', 'GRAND_AVG', thisExp, epIdx, nTotal);
        dataGA = dataGA(:,:,1:epIdx);
        epoch  = epoch(1:epIdx);
        nTotal = epIdx;
    end

    EEGga = eeg_emptyset();
    EEGga.data     = dataGA;
    EEGga.nbchan   = nbchan;
    EEGga.pnts     = pnts;
    EEGga.trials   = nTotal;
    EEGga.srate    = templateEEG.srate;
    EEGga.times    = templateEEG.times;
    EEGga.xmin     = templateEEG.xmin;
    EEGga.xmax     = templateEEG.xmax;
    EEGga.chanlocs = templateEEG.chanlocs;
    EEGga.epoch    = epoch;

    % Store bookkeeping for plotting
    EEGga.etc.grandavg.subjects          = subjIDs;
    EEGga.etc.grandavg.trialCounts_4cond = trialCounts;
    EEGga.etc.grandavg.nSubTotal         = nSubTotal;

    EEGga.etc.grandavg.nTrialLocalStd  = nTrialLocalStd;
    EEGga.etc.grandavg.nTrialLocalDev  = nTrialLocalDev;
    EEGga.etc.grandavg.nTrialGlobalStd = nTrialGlobalStd;
    EEGga.etc.grandavg.nTrialGlobalDev = nTrialGlobalDev;

    EEGga.setname = sprintf('Grand_avg_LG_%s', thisExp);
    EEGga = eeg_checkset(EEGga); % IMPORTANT: no 'eventconsistency' here

    fprintf('   [%s] Nsub=%d | Local trials Std=%d Dev=%d | Global trials Std=%d Dev=%d\n', ...
        thisExp, nSubTotal, nTrialLocalStd, nTrialLocalDev, nTrialGlobalStd, nTrialGlobalDev);

    %% Save grand dataset
    pop_saveset(EEGga, 'filename', [EEGga.setname '.set'], 'filepath', grandFold);

    %% Save evoution plots 
    lg_export_timtopo_subject_and_grand(string(outpath), string(grandFold), string(thisExp), targets, string(skipNow), optsTT, EEGga);

    %% ===================== FIGURES =====================

    % Bekinschtein
    opts = struct();
    opts.topChan     = chans2plot{1};
    opts.botChan     = chans2plot{2};
    opts.winLocal    = [120 220];
    opts.winGlobal   = [300 600];
    opts.tlim        = [-2500 700]; % ms
    opts.ylim        = [-8 8];
    opts.topoGrid    = 80;
    opts.topoShading = 'interp';
    opts.topoContours= 8;
    opts.caxis       = [-1.4 1.4];
    opts.outPng      = fullfile(grandFold, sprintf('Bekinschtein_grand_avg_%s.png', thisExp));
    plot_bekinschtein(EEGga, 'GRAND_AVG', thisExp, opts);

    % Mazancieux
    opts2 = struct();
    opts2.tlim   = [-2500 700];
    opts2.ylim   = [-8 8];
    opts2.outPng = fullfile(grandFold, sprintf('Mazancieux_grand_avg_%s.png', thisExp));
    plot_mz4(EEGga, 'GRAND_AVG', thisExp, opts2);


    % Statistical analyses are performed later in B7_EEG_stats.m.

end % exp loop

disp('Grand averages created + figures saved.');

%% ===== Local helper: reorder dataset channels to template label order =====
function EEGout = reorder_to_template(EEGin, templateLabels)
curLabels = {EEGin.chanlocs.labels};
[~, idx] = ismember(templateLabels, curLabels);
if any(idx == 0)
    missing = templateLabels(idx==0);
    error('Missing channels: %s', strjoin(missing, ', '));
end
EEGout = EEGin;
EEGout.data     = EEGin.data(idx,:,:);
EEGout.chanlocs = EEGin.chanlocs(idx);
EEGout.nbchan   = numel(idx);
EEGout = eeg_checkset(EEGout);
end

%% ===== helper: robust label-to-string-array converter =====
function S = to_strlist(v)
if isstring(v)
    S = v;
elseif ischar(v)
    S = string(v);
elseif iscell(v)
    S = strings(0,1);
    for i=1:numel(v)
        if isstring(v{i}) || ischar(v{i})
            S(end+1,1) = string(v{i}); %#ok<AGROW>
        end
    end
else
    S = strings(0,1);
end
S = strtrim(regexprep(S,'\s+',' '));
end
