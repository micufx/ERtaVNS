clc; clear; close all;

%%
% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% EEG data preparation (Local–Global paradigm)
% This code prepares EEG data for preprocessing. Removes bad channels and
% applies ICA decomposition (on continuous data) for later artifact removal.
%
% Miguel Contreras-Altamirano, 2025

%% Paths & file list
% PILOT SETTING: process one visit at a time. Visit 1 and Visit 2 are kept
% in separate output roots so they can never overwrite each other.
analysisVisit = 1;   % <-- change to 2 when processing Visit 2
PROJECT_ROOT = fileparts(LG_THIS_SCRIPT_DIR);
mainpath = 'L:\Downloads\eeglab2023.0\';  % EEGLAB folder (edit if needed)
path     = fullfile(PROJECT_ROOT,'Experiments','Data');
outpath  = fullfile(PROJECT_ROOT,'Analysis_Output',sprintf('Visit%d',analysisVisit));
if exist(outpath,'dir')~=7, mkdir(outpath); end

% ===================== RUN CONTROL =====================
prepareTasks = true;   % Set false to keep existing E1/E2 ICA and prepare control only.
prepareControl = true; % Also prepare a separate same-visit control ICA set.
startSubIndex = 1;   % <-- Set this to 24 if you want to start from S24 (based on sorted subNames)
% ============================================================

% Add EEGLAB properly
if exist('eeglab','file') ~= 2
    addpath(mainpath);
end

% Current acquisition layout: Experiments/Data/S01/Visit1/S01V1E1.vhdr etc.
% E1+E2 are concatenated ONLY within the selected visit for ICA. Different visits
% are not concatenated because they are separate recording sessions.
% We'll enumerate those folders and, per subject, process E1+E2 *concatenated* to get a single ICA/bad-channel set.

% Event relabeling is NOT needed at the ICA-preparation stage.
% The corrected trigger labels are applied later in B3_EEG_preprocessing_LG.

% ICA filter band (continuous copy used just to learn ICA)
HPF_ica = 1;       % Hz
LPF_ica = 30;      % Hz
thresh_ICA = 5;    % SD threshold for epoch rejection ICA
from = -0.50;      % epoch ICA (sec)
to = 0.50;         % epoch ICA (sec)

% RMS window & thresholds (same plotting logic you used)
sec    = 30;            % seconds per RMS window
thresh_chan = 1;        % SD above mean RMS for bad-channel marking (same as your style) %#ok<NASGU>

%% Outputs we keep (LG-suffixed to avoid confusion with your RP study)
ID_LG           = {};
bad_chans_LG    = {};
No_bad_chans_LG = {};

% Track which subjects were actually finalized in THIS run (so we only overwrite those rows)
processedThisRun = false(1,0);

%% EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% --- Find participant folders S01, S02, ... ---
dSubs  = dir(path);
isSub  = [dSubs.isdir] & ~ismember({dSubs.name},{'.','..'});
allDir = {dSubs(isSub).name};

% Keep only names that look like Sxx (optional, but helpful)
subNames_S = allDir(~cellfun('isempty', regexp(allDir, '^S\d+', 'once')));

% --- Option to include pilot folders too (non-S## names) ---
% resp = questdlg('Also include pilot folders (non S## names)?','Include pilots','No','Yes','No');
% includePilots = strcmpi(resp,'Yes');
includePilots = false;
if includePilots
    pilotNames = setdiff(allDir, subNames_S);
    subNames   = [subNames_S, pilotNames];  % order: S## first, then pilots
else
    subNames   = subNames_S;
end

if isempty(subNames)
    error('No participant folders found under: %s', path);
end

% ===================== NEW: Stable indexing (S01->1, S24->24, etc.) =====================
% Sort so indices are stable and reproducible
subNames = sort(subNames);

% Guard start index
if startSubIndex < 1 || startSubIndex > numel(subNames)
    error('startSubIndex=%d is out of range. Must be between 1 and %d.', startSubIndex, numel(subNames));
end

% Resize run tracker now that subNames length is known
processedThisRun = false(1, numel(subNames));

% If an existing Info_EEG_LG.mat exists, load it so we preserve previous rows AND full table columns
infoMat = fullfile(outpath, 'Info_EEG_LG.mat');
Sfull_prev = struct();         % full MAT contents (so we can preserve all variables on save)
T_LG_prev_full = table();      % full existing table (all columns)
if exist(infoMat, 'file')
    Sfull_prev = load(infoMat);   % <-- load EVERYTHING (important to preserve all vars)

    % --- Keep full previous T_LG (all columns) ---
    if isfield(Sfull_prev,'T_LG') && istable(Sfull_prev.T_LG) && ~isempty(Sfull_prev.T_LG)
        T_LG_prev_full = Sfull_prev.T_LG;
    end

    % --- Prefer existing cell arrays if present ---
    if isfield(Sfull_prev,'ID_LG') && ~isempty(Sfull_prev.ID_LG)
        ID_LG = Sfull_prev.ID_LG;
    end
    if isfield(Sfull_prev,'bad_chans_LG') && ~isempty(Sfull_prev.bad_chans_LG)
        bad_chans_LG = Sfull_prev.bad_chans_LG;
    end
    if isfield(Sfull_prev,'No_bad_chans_LG') && ~isempty(Sfull_prev.No_bad_chans_LG)
        No_bad_chans_LG = Sfull_prev.No_bad_chans_LG;
    end

    % --- Fallback: reconstruct arrays from T_LG if the arrays are missing/empty ---
    % This prevents older rows from being overwritten as [] when resuming from a later subject.
    if istable(T_LG_prev_full) && ~isempty(T_LG_prev_full)
        if isempty(ID_LG) && ismember('Subject_ID', T_LG_prev_full.Properties.VariableNames)
            ID_LG = table_col_to_cell(T_LG_prev_full.Subject_ID);
        end
        if isempty(bad_chans_LG) && ismember('Bad_channels', T_LG_prev_full.Properties.VariableNames)
            bad_chans_LG = table_col_to_cell(T_LG_prev_full.Bad_channels);
        end
        if isempty(No_bad_chans_LG) && ismember('Num_Bad_channels', T_LG_prev_full.Properties.VariableNames)
            No_bad_chans_LG = table_col_to_cell(T_LG_prev_full.Num_Bad_channels);
        end
    end
end

% Normalize to row cell arrays (safer for indexed assignments)
ID_LG           = force_row_cell(ID_LG);
bad_chans_LG    = force_row_cell(bad_chans_LG);
No_bad_chans_LG = force_row_cell(No_bad_chans_LG);

% Pre-allocate to full length so row numbers match subject index
if numel(ID_LG) < numel(subNames)
    ID_LG{numel(subNames)} = [];           %#ok<CCAT>
end
if numel(bad_chans_LG) < numel(subNames)
    bad_chans_LG{numel(subNames)} = [];
end
if numel(No_bad_chans_LG) < numel(subNames)
    No_bad_chans_LG{numel(subNames)} = [];
end

% Fill the Subject_ID column according to the folder ordering
for ii = 1:numel(subNames)
    ID_LG{ii} = subNames{ii};
end
% =======================================================================================

% Process each participant; concatenate E1+E2 (if both exist)
if prepareTasks
for sub = startSubIndex:length(subNames)

    participant = subNames{sub};                                 % e.g., 'S01' or a pilot name
    subFolder   = fullfile(path, participant, sprintf('Visit%d',analysisVisit));
    out_subfold = fullfile(outpath, participant);
    if ~exist(out_subfold,'dir'), mkdir(out_subfold); end

    fprintf('\n[%d/%d] %s\n', sub, numel(subNames), participant);

    % Locate current EyeLink-compatible BrainVision names (allow rerun suffix A/B/...).
    vhdrE1 = dir(fullfile(subFolder, sprintf('%sV%dE1*.vhdr',participant,analysisVisit)));
    vhdrE2 = dir(fullfile(subFolder, sprintf('%sV%dE2*.vhdr',participant,analysisVisit)));

    if isempty(vhdrE1) && isempty(vhdrE2)
        warning('No %sV%dE1/E2 BrainVision files found in %s. Skipping participant.', participant, analysisVisit, subFolder);
        continue;
    end

    % Guard against missing paired .eeg for each .vhdr
    vhdrList = {};  % filenames only
    vhdrE1 = prefer_exact_stem(vhdrE1,sprintf('%sV%dE1.vhdr',participant,analysisVisit));
    vhdrE2 = prefer_exact_stem(vhdrE2,sprintf('%sV%dE2.vhdr',participant,analysisVisit));
    if ~isempty(vhdrE1)
        f1_full = fullfile(subFolder, vhdrE1(1).name);
        if exist(strrep(f1_full,'.vhdr','.eeg'),'file')
            vhdrList{end+1} = vhdrE1(1).name;  % filename only
        else
            warning('Missing binary .eeg for %s (skipping E1 for %s)', vhdrE1(1).name, participant);
        end
    end
    if ~isempty(vhdrE2)
        f2_full = fullfile(subFolder, vhdrE2(1).name);
        if exist(strrep(f2_full,'.vhdr','.eeg'),'file')
            vhdrList{end+1} = vhdrE2(1).name;  % filename only
        else
            warning('Missing binary .eeg for %s (skipping E2 for %s)', vhdrE2(1).name, participant);
        end
    end

    if isempty(vhdrList)
        warning('No loadable E1/E2 BrainVision pair for %s. Skipping.', participant);
        continue;
    end

    %% Load raw continuous BrainVision file(s) and add online reference
    % Load first experiment available, then merge the second (if present) to concatenate E1+E2
    EEG = pop_loadbv(subFolder, vhdrList{1}, [], []);   % <-- path + filename (fixed)

    % Online reference --- After pop_loadbv(...) ---
    origRef = 'FCz';   % <- Used online

    % 1) Set the EEGLAB reference tag (not always filled by importer)
    EEG.ref = origRef;

    % 2) Keep an explicit record in .etc (survives all processing)
    EEG.etc.initial_reference.label  = origRef;
    EEG.etc.initial_reference.kind   = 'hardware';      % or 'software' if virtual
    EEG.etc.initial_reference.source = 'BrainVision online reference';

    % 3) Add a human-readable note in the dataset comments
    EEG.comments = char(EEG.comments, ...
        sprintf('[%s] Initial hardware reference: %s', datestr(now,'yyyy-mm-dd HH:MM'), origRef));
    EEG = eeg_checkset(EEG);

    %% Channel locations
    EEG = pop_chanedit(EEG, 'lookup',[mainpath, 'plugins\dipfit5.4\standard_BESA\standard-10-5-cap385.elp']); % Channel info
    EEG = eeg_checkset(EEG, 'eventconsistency'); % Event consistency


    %% Concatenating EEG data sets (Experiment 1 & 2)

    if numel(vhdrList) == 2

        % load second and merge (concatenate in time)
        EEG2 = pop_loadbv(subFolder, vhdrList{2}, [], []);

        % Online reference --- After pop_loadbv(...) ---
        EEG2.ref = origRef;  % 1) Set the EEGLAB reference tag (not always filled by importer)
        EEG2.etc.initial_reference.label  = origRef;  % 2) Keep an explicit record in .etc (survives all processing)
        EEG2.etc.initial_reference.kind   = 'hardware';      % or 'software' if virtual
        EEG2.etc.initial_reference.source = 'BrainVision online reference';
        EEG2.comments = char(EEG2.comments, ...  % 3) Add a human-readable note in the dataset comments
            sprintf('[%s] Initial hardware reference: %s', datestr(now,'yyyy-mm-dd HH:MM'), origRef));
        EEG2 = eeg_checkset(EEG2);

        % Channel locations
        EEG2 = pop_chanedit(EEG2, 'lookup',[mainpath, 'plugins\dipfit5.4\standard_BESA\standard-10-5-cap385.elp']); % Channel info
        EEG2 = eeg_checkset(EEG2, 'eventconsistency'); % Event consistency

        % Store and merge via EEGLAB (concatenate datasets)
        [ALLEEG, EEG]  = eeg_store(ALLEEG, EEG,  1);
        [ALLEEG, EEG2] = eeg_store(ALLEEG, EEG2, 2);
        EEG = pop_mergeset(ALLEEG, [1 2], 0);
        EEG = eeg_checkset(EEG);
        ALLEEG = []; CURRENTSET = 1;  %#ok<NASGU> % clean store to keep memory simple

    end
    eeglab redraw % Updating GUI

    % Build a label for plots/files: participant + (E1+E2 or single)
    if numel(vhdrList) == 2
        baseName = [participant, '_E1E2'];
    else
        [~, stemSingle] = fileparts(vhdrList{1});
        baseName = stemSingle;
    end

    %% Power-line cleanup (after loading/merging & chanlocs)
    % Notch filter(multitaper regression) → minimal waveform distortion for
    % ERP / Pre-post line-noise PSD report
    [EEG, figPSD] = clean_line_noise(EEG, ...
    'linefreq',50, ...                % 50 EU / 60 US
    'useCleanline',true, ...          % try CleanLine first
    'outPng', fullfile(out_subfold, ['PSD_LineNoise_', baseName, '.png'])); %#ok<NASGU>

    % Optional
    %  Harmonics up to Nyquist: maxHarmonics = floor((EEG.srate/2)/50)
    %  (e.g., at 500 Hz sampling: up to 200 Hz → 4 harmonics).

    % 'maxHarmonics',4, ...             % targets 50,100,150,200 (if below Nyquist)
    % 'bandwidthHz',3, ...              % slightly wider than default 2 Hz

    %% Bad-channel screening via RMS
    nchan_ori = EEG.nbchan;
    rms_all   = std(double(EEG.data)');  % RMS per channel (across time)

    % MAD = median absolute deviation method
    USE_ROBUST = false;    % <— set false to keep mean+SD style

    if USE_ROBUST
        medR = median(rms_all);
        madR = median(abs(rms_all - medR));
        robustSD = 1.4826 * madR;                 % ~= SD if normal
        Z = (rms_all - medR) / robustSD;
        zCut = 3.0;                               % typical robust cutoff
        ind_bad = find(Z > zCut);
    else

        % Mean + SD method (default)
        mean_RMS  = mean(rms_all);
        sd_RMS  = std(rms_all);
        num_SD = 1.5;                                 % 1.0 SD (current), or 2.0 more conservative
        thresh_chan = mean_RMS + num_SD*sd_RMS;       % Define threshold for bad channel marking
        ind_bad   = find(rms_all > thresh_chan);
    end


    % RMS over time (30 s windows)
    LeWin     = EEG.srate*sec;
    lastIdx   = size(EEG.data,2)-LeWin+1;
    if lastIdx < 1
        warning('Recording shorter than %d s window; skipping RMS-over-time for %s', sec, baseName);
        rms_t = zeros(EEG.nbchan,0);
    else
        idx_loop  = 1:LeWin:lastIdx;
        rms_t     = zeros(EEG.nbchan, numel(idx_loop));
        row_count = 1;
        for idx = idx_loop
            seg = double(EEG.data(:, idx:(idx+LeWin-1)));
            rms_t(:, row_count) = std(seg, [], 2);
            row_count = row_count + 1;
        end
    end

    % Figure
    fig_bad_chan = figure('units','normalized','outerposition',[0 0 1 1]);

    subtightplot(2, 2, [1, 3]);
    topoplot(rms_all, EEG.chanlocs, 'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo, ...
        'whitebk','on','shading','interp');
    title('RMS of Electrodes Over Time');

    subplot(2,2,4);
    plot(rms_all,'k','LineWidth',1); hold on;
    plot(repmat(thresh_chan,1,EEG.nbchan),'r--','LineWidth',2);
    title(['Outliers: ', mat2str(ind_bad)], 'Interpreter','none');
    subtitle('Threshold: [Mean RMS + 1.5 SD]');
    xlabel('Electrodes'); ylabel('RMS (\muV)'); grid on; axis tight;

    ax1 = subplot(2,2,2);
    if ~isempty(rms_t)
        imagesc(rms_t); colorbar; colormap(ax1, turbo);
        xlabel('Time (s)'); ylabel('Electrodes');
        title('Root Mean Square [RMS]');
        subtitle(['Windows of ', num2str(sec), ' s']);
    else
        text(0.5,0.5,'Recording too short for 30s RMS windows','HorizontalAlignment','center');
        axis off;
    end

    sgtitle(['Electrode artefacts: [', baseName, ']'], 'Interpreter','none');

    % Mark/remove
    EEG.badchans = sort(ind_bad);
    no_bad_chans = numel(EEG.badchans);
    if no_bad_chans>0
        EEG = pop_select(EEG, 'rmchannel', EEG.badchans);
        EEG = eeg_checkset(EEG);
    end
    disp(['Original #ch: ', num2str(nchan_ori), '. Bad channels: ', num2str(no_bad_chans)]);

    % Save bad-channel plot
    saveas(fig_bad_chan, fullfile(out_subfold, ['Bad_chan_', baseName, '.jpg']));
    close(fig_bad_chan);

    %% ---- Learn ICA on continuous, 1–40 Hz (no epoching; cleaner for auditory) ----
    %% === OPTION A: Learn ICA on cotinuous data (toggle this block on to compare) ===
    EEG_ica = EEG;                                     % start from continuous
    EEG_ica = pop_eegfiltnew(EEG_ica,'locutoff',HPF_ica,'hicutoff',LPF_ica);
    EEG_ica = eeg_checkset(EEG_ica); % check dataset
    EEG_ica = pop_runica(EEG_ica, 'icatype','runica','extended',1,'interrupt','on');
    EEG_ica = eeg_checkset(EEG_ica); % check dataset

    %% === OPTION B: Learn ICA on 1-s pseudo-epochs (toggle this block on to compare) ===
    % EEG_ica = EEG;                                     % start from continuous
    % EEG_ica = pop_eegfiltnew(EEG_ica,'locutoff',HPF_ica,'hicutoff',LPF_ica);
    %
    % EEG_ica = eeg_regepochs(EEG_ica, 'recurrence',1, 'limits',[0 1], 'rmbase',NaN); % reject 1sec consecutive epochs
    %
    % EEG_ica = eeg_checkset(EEG_ica); % check dataset
    % EEG_ica = pop_runica(EEG_ica,'icatype','runica','extended',1,'interrupt','off');
    % EEG_ica = eeg_checkset(EEG_ica); % check dataset

    %% Store weights back to the (same-channel) continuous EEG
    EEG.icaact      = [];                      % keep memory small / or EEG.icaact = EEG_ica.icaact;
    EEG.icawinv     = EEG_ica.icawinv;
    EEG.icasphere   = EEG_ica.icasphere;
    EEG.icaweights  = EEG_ica.icaweights;
    EEG.icachansind = EEG_ica.icachansind;
    EEG.setname     = [EEG.setname, '_ica'];

    clear EEG_ica;
    eeglab redraw % Updating GUI

    % --- Save ONE dataset per participant (concatenated); name by participant ---
    pop_saveset(EEG, 'filename', ['ica_', participant, '.set'], 'filepath', out_subfold);

    % Bookkeeping for the LG table (ONE row per participant)
    % ===================== NEW: write into the correct ROW =====================
    ID_LG{sub}           = participant;
    bad_chans_LG{sub}    = EEG.badchans;
    No_bad_chans_LG{sub} = no_bad_chans;
    processedThisRun(sub) = true;   % mark as updated in this run
    % ==========================================================================

    disp([participant, ' finalized!']);

end

end % prepareTasks

%% Table info (LG)
% ===================== NEW: preserve ALL existing T_LG columns =====================
% Start from previous full table if it exists; otherwise create a new one
if exist('T_LG_prev_full','var') && istable(T_LG_prev_full) && ~isempty(T_LG_prev_full)
    T_LG = T_LG_prev_full;
else
    T_LG = table();
end

% Ensure Subject_ID column exists
if ~ismember('Subject_ID', T_LG.Properties.VariableNames)
    T_LG.Subject_ID = cell(height(T_LG),1);
end

% Ensure the 3 columns we manage exist (create if missing)
if ~ismember('Bad_channels', T_LG.Properties.VariableNames)
    T_LG.Bad_channels = repmat({[]}, height(T_LG), 1);
end
if ~ismember('Num_Bad_channels', T_LG.Properties.VariableNames)
    T_LG.Num_Bad_channels = repmat({[]}, height(T_LG), 1);
end

% Update ONLY subjects processed in this run (do not touch others)
for sub = 1:numel(subNames)
    if ~processedThisRun(sub)
        continue;
    end

    participant = subNames{sub};

    % Find row by Subject_ID (robust)
    row = [];
    if height(T_LG) > 0
        try
            row = find(strcmpi(string(T_LG.Subject_ID), string(participant)), 1);
        catch
            % If Subject_ID has mixed types, normalize through helper
            subjCol = table_col_to_cell(T_LG.Subject_ID);
            row = find(strcmpi(string(subjCol), string(participant)), 1);
        end
    end

    % If row doesn't exist, append a new empty row preserving all columns
    if isempty(row)
        T_LG = append_empty_row_like(T_LG);
        row = height(T_LG);
    end

    % Write/overwrite ONLY the fields handled by this script
    T_LG.Subject_ID{row} = participant;
    T_LG.Bad_channels{row} = bad_chans_LG{sub};
    T_LG.Num_Bad_channels{row} = No_bad_chans_LG{sub};
end

% Keep the lightweight arrays synchronized (optional but useful for compatibility)
% Build arrays from the first numel(subNames) rows by subject folder order when possible
ID_LG_out           = ID_LG;
bad_chans_LG_out    = bad_chans_LG;
No_bad_chans_LG_out = No_bad_chans_LG;

ID_LG           = ID_LG_out;
bad_chans_LG    = bad_chans_LG_out;
No_bad_chans_LG = No_bad_chans_LG_out;
% ================================================================================

if ~exist(outpath,'dir'), mkdir(outpath); end

% Preserve all other variables already inside Info_EEG_LG.mat
Ssave = struct();
if exist('Sfull_prev','var') && isstruct(Sfull_prev) && ~isempty(fieldnames(Sfull_prev))
    Ssave = Sfull_prev;
end
Ssave.T_LG = T_LG;
Ssave.ID_LG = ID_LG;
Ssave.bad_chans_LG = bad_chans_LG;
Ssave.No_bad_chans_LG = No_bad_chans_LG;

save(fullfile(outpath, 'Info_EEG_LG.mat'), '-struct', 'Ssave');
disp('Saved Info_EEG_LG.mat');

% Write Excel export
excelFile = fullfile(outpath,'Info_EEG_LG.xlsx');
write_TLG_excel_exact(T_LG, excelFile);
fprintf('Excel export written: %s\n', excelFile);

%% Control ICA preparation (same settings; separate decomposition)
% Preserve the established E1/E2 ICA and its Info_EEG_LG table entries.
if prepareControl
    for sub = startSubIndex:numel(subNames)
        participant = subNames{sub};
        subFolder = fullfile(path,participant,sprintf('Visit%d',analysisVisit));
        out_subfold = fullfile(outpath,participant);
        d = dir(fullfile(subFolder,sprintf('%sV%dC*.vhdr',participant,analysisVisit)));
        exact = strcmpi({d.name},sprintf('%sV%dC.vhdr',participant,analysisVisit));
        if any(exact), d=d(exact); end
        if isempty(d), continue; end
        try
            assert(numel(d)==1,'Multiple control recordings: keep one exact SxxV#C stem.');
            assert(isfile(fullfile(subFolder,strrep(d.name,'.vhdr','.eeg'))),'Control EEG binary is missing.');
            if ~isfolder(out_subfold), mkdir(out_subfold); end
            prepare_control_ica(subFolder,d.name,out_subfold,participant,mainpath,HPF_ica,LPF_ica,sec);
        catch ME
            warning('B2:ControlFailed','%s V%d control preparation failed: %s',participant,analysisVisit,ME.message);
        end
    end
end

%% ========================= LOCAL HELPERS =========================
function c = table_col_to_cell(x)
% Convert a table column (cell/string/numeric/etc.) to a row cell array
if iscell(x)
    c = x(:)';
elseif isstring(x)
    c = cellstr(x(:))';
elseif isnumeric(x) || islogical(x)
    c = num2cell(x(:))';
else
    try
        c = num2cell(x(:))';
    catch
        c = {};
    end
end
end

function c = force_row_cell(c)
% Ensure row cell array for safe indexed assignment with {idx}
if isempty(c)
    c = {};
    return;
end
if ~iscell(c)
    c = table_col_to_cell(c);
else
    c = c(:)';
end
end

function T = append_empty_row_like(T)
% Append one empty row preserving all existing columns and variable types
if ~istable(T)
    T = table();
    return;
end

if height(T) == 0
    % Create one-row empties for each variable based on current variable classes
    names = T.Properties.VariableNames;
    for k = 1:numel(names)
        v = T.(names{k});
        T.(names{k})(1,1) = make_empty_scalar_like(v);
    end
else
    % Generic and robust for most table variable types
    T = [T; T(1,:)]; %#ok<AGROW>
    r = height(T);
    names = T.Properties.VariableNames;
    for k = 1:numel(names)
        col = T.(names{k});
        try
            if iscell(col)
                T.(names{k}){r} = [];
            elseif isstring(col)
                T.(names{k})(r) = "";
            elseif isnumeric(col)
                T.(names{k})(r) = NaN;
            elseif islogical(col)
                T.(names{k})(r) = false;
            elseif ischar(col)
                % char matrix columns are uncommon in tables; leave as copied row
            else
                % leave as copied row if type-specific assignment is unclear
            end
        catch
            % leave copied value if assignment fails
        end
    end
end
end

function x = make_empty_scalar_like(v)
% Helper used if table has 0 rows (rare case)
if iscell(v)
    x = {[]};
elseif isstring(v)
    x = "";
elseif isnumeric(v)
    x = NaN;
elseif islogical(v)
    x = false;
else
    try
        x = v([]);
    catch
        x = [];
    end
end
end

function d = prefer_exact_stem(d, exactName)
if isempty(d), return; end
idx=find(strcmpi({d.name},exactName),1);
if ~isempty(idx), d=d([idx setdiff(1:numel(d),idx,'stable')]); end
end

function prepare_control_ica(subFolder,fileName,out_subfold,participant,mainpath,HPF_ica,LPF_ica,sec)
% Same continuous line-noise cleanup, RMS screen and ICA fit as task preparation.
vhdrList={fileName};
    %% Load raw continuous BrainVision file(s) and add online reference
    % Load first experiment available, then merge the second (if present) to concatenate E1+E2
    EEG = pop_loadbv(subFolder, vhdrList{1}, [], []);   % <-- path + filename (fixed)

    % Online reference --- After pop_loadbv(...) ---
    origRef = 'FCz';   % <- Used online

    % 1) Set the EEGLAB reference tag (not always filled by importer)
    EEG.ref = origRef;

    % 2) Keep an explicit record in .etc (survives all processing)
    EEG.etc.initial_reference.label  = origRef;
    EEG.etc.initial_reference.kind   = 'hardware';      % or 'software' if virtual
    EEG.etc.initial_reference.source = 'BrainVision online reference';

    % 3) Add a human-readable note in the dataset comments
    EEG.comments = char(EEG.comments, ...
        sprintf('[%s] Initial hardware reference: %s', datestr(now,'yyyy-mm-dd HH:MM'), origRef));
    EEG = eeg_checkset(EEG);

    %% Channel locations
    EEG = pop_chanedit(EEG, 'lookup',[mainpath, 'plugins\dipfit5.4\standard_BESA\standard-10-5-cap385.elp']); % Channel info
    EEG = eeg_checkset(EEG, 'eventconsistency'); % Event consistency


    [~,baseName]=fileparts(fileName);

    % Keep control preparation separate from the working E1/E2 decomposition.
    originalLabels = {EEG.chanlocs.labels};
    %% Power-line cleanup (after loading/merging & chanlocs)
    % Notch filter(multitaper regression) → minimal waveform distortion for
    % ERP / Pre-post line-noise PSD report
    [EEG, figPSD] = clean_line_noise(EEG, ...
    'linefreq',50, ...                % 50 EU / 60 US
    'useCleanline',true, ...          % try CleanLine first
    'outPng', fullfile(out_subfold, ['PSD_LineNoise_', baseName, '.png'])); %#ok<NASGU>

    % Optional
    %  Harmonics up to Nyquist: maxHarmonics = floor((EEG.srate/2)/50)
    %  (e.g., at 500 Hz sampling: up to 200 Hz → 4 harmonics).

    % 'maxHarmonics',4, ...             % targets 50,100,150,200 (if below Nyquist)
    % 'bandwidthHz',3, ...              % slightly wider than default 2 Hz

    %% Bad-channel screening via RMS
    nchan_ori = EEG.nbchan;
    rms_all   = std(double(EEG.data)');  % RMS per channel (across time)

    % MAD = median absolute deviation method
    USE_ROBUST = false;    % <— set false to keep mean+SD style

    if USE_ROBUST
        medR = median(rms_all);
        madR = median(abs(rms_all - medR));
        robustSD = 1.4826 * madR;                 % ~= SD if normal
        Z = (rms_all - medR) / robustSD;
        zCut = 3.0;                               % typical robust cutoff
        ind_bad = find(Z > zCut);
    else

        % Mean + SD method (default)
        mean_RMS  = mean(rms_all);
        sd_RMS  = std(rms_all);
        num_SD = 1.5;                                 % 1.0 SD (current), or 2.0 more conservative
        thresh_chan = mean_RMS + num_SD*sd_RMS;       % Define threshold for bad channel marking
        ind_bad   = find(rms_all > thresh_chan);
    end


    % RMS over time (30 s windows)
    LeWin     = EEG.srate*sec;
    lastIdx   = size(EEG.data,2)-LeWin+1;
    if lastIdx < 1
        warning('Recording shorter than %d s window; skipping RMS-over-time for %s', sec, baseName);
        rms_t = zeros(EEG.nbchan,0);
    else
        idx_loop  = 1:LeWin:lastIdx;
        rms_t     = zeros(EEG.nbchan, numel(idx_loop));
        row_count = 1;
        for idx = idx_loop
            seg = double(EEG.data(:, idx:(idx+LeWin-1)));
            rms_t(:, row_count) = std(seg, [], 2);
            row_count = row_count + 1;
        end
    end

    % Figure
    fig_bad_chan = figure('units','normalized','outerposition',[0 0 1 1]);

    subtightplot(2, 2, [1, 3]);
    topoplot(rms_all, EEG.chanlocs, 'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo, ...
        'whitebk','on','shading','interp');
    title('RMS of Electrodes Over Time');

    subplot(2,2,4);
    plot(rms_all,'k','LineWidth',1); hold on;
    plot(repmat(thresh_chan,1,EEG.nbchan),'r--','LineWidth',2);
    title(['Outliers: ', mat2str(ind_bad)], 'Interpreter','none');
    subtitle('Threshold: [Mean RMS + 1.5 SD]');
    xlabel('Electrodes'); ylabel('RMS (\muV)'); grid on; axis tight;

    ax1 = subplot(2,2,2);
    if ~isempty(rms_t)
        imagesc(rms_t); colorbar; colormap(ax1, turbo);
        xlabel('Time (s)'); ylabel('Electrodes');
        title('Root Mean Square [RMS]');
        subtitle(['Windows of ', num2str(sec), ' s']);
    else
        text(0.5,0.5,'Recording too short for 30s RMS windows','HorizontalAlignment','center');
        axis off;
    end

    sgtitle(['Electrode artefacts: [', baseName, ']'], 'Interpreter','none');

    % Mark/remove
    EEG.badchans = sort(ind_bad);
    no_bad_chans = numel(EEG.badchans);
    if no_bad_chans>0
        EEG = pop_select(EEG, 'rmchannel', EEG.badchans);
        EEG = eeg_checkset(EEG);
    end
    disp(['Original #ch: ', num2str(nchan_ori), '. Bad channels: ', num2str(no_bad_chans)]);

    % Save bad-channel plot
    saveas(fig_bad_chan, fullfile(out_subfold, ['Bad_chan_', baseName, '.jpg']));
    close(fig_bad_chan);

    %% ---- Learn ICA on continuous, 1–40 Hz (no epoching; cleaner for auditory) ----
    %% === OPTION A: Learn ICA on cotinuous data (toggle this block on to compare) ===
    EEG_ica = EEG;                                     % start from continuous
    EEG_ica = pop_eegfiltnew(EEG_ica,'locutoff',HPF_ica,'hicutoff',LPF_ica);
    EEG_ica = eeg_checkset(EEG_ica); % check dataset
    EEG_ica = pop_runica(EEG_ica, 'icatype','runica','extended',1,'interrupt','on');
    EEG_ica = eeg_checkset(EEG_ica); % check dataset

    %% === OPTION B: Learn ICA on 1-s pseudo-epochs (toggle this block on to compare) ===
    % EEG_ica = EEG;                                     % start from continuous
    % EEG_ica = pop_eegfiltnew(EEG_ica,'locutoff',HPF_ica,'hicutoff',LPF_ica);
    %
    % EEG_ica = eeg_regepochs(EEG_ica, 'recurrence',1, 'limits',[0 1], 'rmbase',NaN); % reject 1sec consecutive epochs
    %
    % EEG_ica = eeg_checkset(EEG_ica); % check dataset
    % EEG_ica = pop_runica(EEG_ica,'icatype','runica','extended',1,'interrupt','off');
    % EEG_ica = eeg_checkset(EEG_ica); % check dataset

    %% Store weights back to the (same-channel) continuous EEG
    EEG.icaact      = [];                      % keep memory small / or EEG.icaact = EEG_ica.icaact;
    EEG.icawinv     = EEG_ica.icawinv;
    EEG.icasphere   = EEG_ica.icasphere;
    EEG.icaweights  = EEG_ica.icaweights;
    EEG.icachansind = EEG_ica.icachansind;
    EEG.setname     = [EEG.setname, '_ica'];

    clear EEG_ica;
    eeglab redraw % Updating GUI

    EEG.etc.control_preparation = struct('source_vhdr',vhdrList{1}, ...
        'bad_channel_labels',{originalLabels(ind_bad)},'bad_channel_indices',ind_bad, ...
        'ica_highpass_hz',HPF_ica,'ica_lowpass_hz',LPF_ica, ...
        'original_channel_labels',{originalLabels});
    % --- Save ONE dataset per participant (concatenated); name by participant ---
    pop_saveset(EEG, 'filename', ['ica_', participant, '_C.set'], 'filepath', out_subfold);


end
