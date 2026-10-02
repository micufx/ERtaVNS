clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% GRAND-AVERAGE PUPIL ANALYSIS (Across Participants) — E1 & E2
%
% PURPOSE
%   This script computes grand-average pupil responses across participants
%   using the preprocessed outputs created by the pupil preprocessing pipeline
%   (i.e., the *_E1_pupil_preproc.mat and *_E2_pupil_preproc.mat files).
%
% INPUTS (per participant, per experiment)
%   Expected file name pattern inside each participant folder:
%       <SUBJECT>_<E1/E2>_pupil_preproc.mat
%
%   Required variables inside each MAT file:
%       - dataTableClean : table of retained trials after preprocessing
%           Required columns:
%             * sound   : "Local Standard" or "Local Deviant"
%             * type    : one of:
%                 "Frequent Local Standard"
%                 "Rare Local Standard"
%                 "Frequent Local Deviant"
%                 "Rare Local Deviant"
%
% Local Standard (LS) = FLS + RLS
% Local Deviant (LD) = FLD + RLD
% Global Standard (GS) = FLS + FLD (all frequent patterns)
% Global Deviant (GD) = RLS + RLD (all rare patterns)
%
%             * PUPIL_Z  : cell array, each cell is a 1 x nTime vector
%                          (z-normalized pupil timecourse per trial)
%       - relTimes : 1 x nTime vector of epoch timepoints in ms
%
% GRAND AVERAGING STRATEGY (IMPORTANT)
%   Equal-weight participant grand average:
%     1) Subject mean per condition (mean over trials)
%     2) Grand mean = mean over subject means (each subject contributes equally)
%
%   SEM is computed across participants at each time point (time-resolved SEM).
%
% CONDITIONS / OUTPUT PLOTS
%   1) SOUND plot (Local Standard vs Local Deviant)
%   2) TYPE plot  (Frequent/Rare x Local Std/Dev)
%   3) Bekinschtein-like 2-panel plot in ONE FIGURE:
%        Top   = Local effect (Local Standard vs Local Deviant)
%        Bottom= Global effect (Global Standard vs Global Deviant)
%      Titles include Nsub and total Ntrials per condition.
%
% ADDITIONAL OUTPUT (requested)
%   Stores per-subject pupil trial counts into T_LG (Info_EEG_LG.mat)
%   and writes an Excel export:
%     - Info_EEG_LG.xlsx   (using write_TLG_excel_exact)
%
% AUTHOR / PROJECT
%   Miguel Contreras Altamirano — Pupil preprocessing & grand-average pipeline
% ========================================================================


%% ===================== USER SETTINGS =====================
outpath  = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';
useField = 'PUPIL_Z';   % from your preprocessing outputs

saveFigures  = true;
saveGrandMat = true;

grandOutFolder = fullfile(outpath, 'GRAND_AVG');
if ~exist(grandOutFolder, 'dir'); mkdir(grandOutFolder); end

% ===================== NEW: OPTIONAL EXCLUSION LIST =====================
% Put here any participants you want to exclude from the group analysis.
% Use the folder names exactly as they appear (e.g., "S01", "S02", ...).
excludeSubs = ["S21","S12","S23","S02","S29"];  % example; set to [] if none  excludeSubs = ["S03","S11","S19"];
% ======================================================================= 

%% ===================== LOAD T_LG (for storing summary) =====================
infoMat = fullfile(outpath,'Info_EEG_LG.mat');
if ~exist(infoMat,'file')
    error('Info_EEG_LG.mat not found at: %s', infoMat);
end
load(infoMat,'T_LG');

% Ensure columns exist and are placed after BaselineAnchor_E2_Pupil
T_LG = ensure_pupil_grandavg_columns(T_LG);

%% ===================== FIND PARTICIPANTS =====================
% Robust: look specifically for folders starting with "S" under outpath
if ~exist(outpath,'dir')
    error('outpath does not exist or is not mounted: %s', outpath);
end

dSubs = dir(fullfile(outpath, 'S*'));   % only folders starting with S
dSubs = dSubs([dSubs.isdir]);           % keep only directories

% Remove "." ".." and GRAND_AVG if it matches
names = string({dSubs.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));

% Keep things like S01, S02, S10, etc.
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+', 'once'));
subNames = cellstr(names(keep));

if isempty(subNames)
    error('No participant folders found under: %s  (expected folders like S01, S02, ...)', outpath);
end

% ===================== NEW: APPLY EXCLUSION LIST =====================
if ~isempty(excludeSubs)
    beforeN = numel(subNames);
    subNames = subNames(~ismember(string(subNames), excludeSubs));
    afterN  = numel(subNames);

    fprintf('Excluded %d participant(s): %s\n', beforeN-afterN, strjoin(excludeSubs, ", "));
end
fprintf('Participants included for analysis: %d\n', numel(subNames));

if isempty(subNames)
    error('All participants were excluded. Nothing to analyze.');
end
% =====================================================================

%% ===================== CONDITIONS (MUST MATCH YOUR LABELS) =====================
soundOrder = ["Local Standard", "Local Deviant"];

typeOrder  = [ ...
    "Frequent Local Standard", ...
    "Rare Local Standard", ...
    "Frequent Local Deviant", ...
    "Rare Local Deviant"];

% Derived global effect labels (computed from TYPE if no explicit column exists)
effectOrder = ["Global Standard", "Global Deviant"];

%% ===================== MAIN LOOP: E1 & E2 =====================
GRAND = struct();

for eidx = 1:2
    thisExp = sprintf('E%d', eidx);
    fprintf('\n====================================================\n');
    fprintf('GRAND AVERAGE: %s\n', thisExp);
    fprintf('====================================================\n');

    % Subject-level means (equal-weight grand average)
    subjMeans_sound  = cell(numel(soundOrder), 1);  % each: [nSub x nTime]
    subjMeans_type   = cell(numel(typeOrder), 1);
    subjMeans_effect = cell(numel(effectOrder), 1);

    % For titles (total trials across loaded subjects)
    totalTrials_sound  = zeros(numel(soundOrder), 1);   % [LS; LD]
    totalTrials_effect = zeros(numel(effectOrder), 1);  % [GS; GD]

    relTimes_ref     = [];
    subjCountLoaded  = 0;

    for s = 1:numel(subNames)
        participant = string(subNames{s});
        subFolder   = fullfile(outpath, char(participant));

        matFile = fullfile(subFolder, sprintf('%s_%s_pupil_preproc.mat', participant, thisExp));
        if ~exist(matFile, 'file')
            fprintf('  [%s | %s] MAT not found -> skipping\n', participant, thisExp);
            continue;
        end

        fprintf('\n  [%s | %s] Loading: %s\n', participant, thisExp, matFile);
        load(matFile, 'dataTableClean', 'relTimes');

        if ~exist('dataTableClean','var') || isempty(dataTableClean) || height(dataTableClean)==0
            fprintf('  [%s | %s] dataTableClean empty -> skipping\n', participant, thisExp);
            clear dataTableClean relTimes
            continue;
        end

        % Check time vector consistency
        if isempty(relTimes_ref)
            relTimes_ref = relTimes(:)';
            fprintf('  [%s | %s] relTimes set as reference (%d samples)\n', participant, thisExp, numel(relTimes_ref));
        else
            if numel(relTimes) ~= numel(relTimes_ref) || any(relTimes(:)' ~= relTimes_ref)
                error('Time vector mismatch in %s (%s). Ensure all files use same relTimes.', participant, thisExp);
            end
        end

        subjCountLoaded = subjCountLoaded + 1;
        nRetained = height(dataTableClean);
        fprintf('  [%s | %s] Retained trials in dataTableClean: %d\n', participant, thisExp, nRetained);

        %% ---------- Trial counts for LOCAL (SOUND) ----------
        soundVec = string(dataTableClean.sound);
        nLS = sum(soundVec == soundOrder(1)); % Local Standard
        nLD = sum(soundVec == soundOrder(2)); % Local Deviant
        fprintf('  [%s | %s] LOCAL trial counts -> LS: %d | LD: %d\n', participant, thisExp, nLS, nLD);

        totalTrials_sound(1) = totalTrials_sound(1) + nLS;
        totalTrials_sound(2) = totalTrials_sound(2) + nLD;

        %% ---------- Trial counts for TYPE (sanity check) ----------
        typeVec = string(dataTableClean.type);
        nFLS = sum(typeVec == typeOrder(1));
        nRLS = sum(typeVec == typeOrder(2));
        nFLD = sum(typeVec == typeOrder(3));
        nRLD = sum(typeVec == typeOrder(4));

        fprintf('  [%s | %s] TYPE trial counts -> %s: %d | %s: %d | %s: %d | %s: %d\n', ...
            participant, thisExp, typeOrder(1), nFLS, typeOrder(2), nRLS, typeOrder(3), nFLD, typeOrder(4), nRLD);

        %% ---------- Derive GLOBAL effect labels from TYPE ----------
        % Frequent* = Global Standard; Rare* = Global Deviant
        effectVec = strings(height(dataTableClean),1);
        effectVec(startsWith(typeVec, "Frequent")) = "Global Standard";
        effectVec(startsWith(typeVec, "Rare"))     = "Global Deviant";

        nGS = sum(effectVec == effectOrder(1)); % Global Standard (Frequent)
        nGD = sum(effectVec == effectOrder(2)); % Global Deviant  (Rare)
        fprintf('  [%s | %s] GLOBAL trial counts -> GS: %d | GD: %d\n', participant, thisExp, nGS, nGD);

        totalTrials_effect(1) = totalTrials_effect(1) + nGS;
        totalTrials_effect(2) = totalTrials_effect(2) + nGD;

        %% ---------- Store subject counts into T_LG (requested) ----------
        totalDesign = 270;
        nBad = totalDesign - nRetained;

        rowIdx = find(strcmpi(string(T_LG.Subject_ID), participant), 1);
        if isempty(rowIdx)
            warning('T_LG.Subject_ID does not contain %s -> not storing pupil summary for this subject.', participant);
        else
            % Good/Bad
            T_LG.("Num_good_trials_"+thisExp+"_Pupil")(rowIdx) = nRetained;
            T_LG.("Num_bad_trials_"+thisExp+"_Pupil")(rowIdx)  = nBad;

            % Bekinschtein-style collapsed counts
            T_LG.("nLS_"+thisExp+"_Pupil")(rowIdx) = nLS;
            T_LG.("nLD_"+thisExp+"_Pupil")(rowIdx) = nLD;
            T_LG.("nGS_"+thisExp+"_Pupil")(rowIdx) = nGS;
            T_LG.("nGD_"+thisExp+"_Pupil")(rowIdx) = nGD;

            % Keep the older SOUND counts too (explicitly named as SOUND)
            T_LG.("nStd_sound_"+thisExp+"_Pupil")(rowIdx) = nLS;
            T_LG.("nDev_sound_"+thisExp+"_Pupil")(rowIdx) = nLD;

            % Keep the 4 TYPE counts as well
            T_LG.("nFLS_"+thisExp+"_Pupil")(rowIdx) = nFLS;
            T_LG.("nRLS_"+thisExp+"_Pupil")(rowIdx) = nRLS;
            T_LG.("nFLD_"+thisExp+"_Pupil")(rowIdx) = nFLD;
            T_LG.("nRLD_"+thisExp+"_Pupil")(rowIdx) = nRLD;
        end

        %% ---------- Compute subject means per condition ----------
        % SOUND means per subject
        for c = 1:numel(soundOrder)
            cond = soundOrder(c);
            m = subject_mean_curve(dataTableClean, useField, 'sound', cond);

            if all(isnan(m))
                fprintf('  [%s | %s]   -> No trials for LOCAL "%s" (subject mean = NaNs)\n', participant, thisExp, cond);
            else
                fprintf('  [%s | %s]   -> Computed subject mean for LOCAL "%s"\n', participant, thisExp, cond);
            end
            subjMeans_sound{c}(end+1, :) = m; %#ok<SAGROW>
        end

        % TYPE means per subject
        for c = 1:numel(typeOrder)
            cond = typeOrder(c);
            m = subject_mean_curve(dataTableClean, useField, 'type', cond);

            if all(isnan(m))
                fprintf('  [%s | %s]   -> No trials for TYPE "%s" (subject mean = NaNs)\n', participant, thisExp, cond);
            else
                fprintf('  [%s | %s]   -> Computed subject mean for TYPE "%s"\n', participant, thisExp, cond);
            end
            subjMeans_type{c}(end+1, :) = m; %#ok<SAGROW>
        end

        % GLOBAL effect means per subject (using derived effectVec)
        Ttmp = dataTableClean;
        Ttmp.effect = effectVec;

        for c = 1:numel(effectOrder)
            cond = effectOrder(c);
            m = subject_mean_curve(Ttmp, useField, 'effect', cond);

            if all(isnan(m))
                fprintf('  [%s | %s]   -> No trials for GLOBAL "%s" (subject mean = NaNs)\n', participant, thisExp, cond);
            else
                fprintf('  [%s | %s]   -> Computed subject mean for GLOBAL "%s"\n', participant, thisExp, cond);
            end
            subjMeans_effect{c}(end+1, :) = m; %#ok<SAGROW>
        end

        clear dataTableClean relTimes Ttmp effectVec
    end

    if subjCountLoaded == 0
        warning('No subjects loaded for %s. Skipping.', thisExp);
        continue;
    end

    %% ---------- Console summary ----------
    fprintf('\n==================== SUMMARY %s ====================\n', thisExp);
    fprintf('Subjects loaded (non-empty dataTableClean): %d\n', subjCountLoaded);
    fprintf('Total retained trials -> LS: %d | LD: %d | GS: %d | GD: %d\n', ...
        totalTrials_sound(1), totalTrials_sound(2), totalTrials_effect(1), totalTrials_effect(2));
    fprintf('====================================================\n');

    %% ---------- Grand mean + SEM across participants ----------
    fprintf('[%s] Computing grand mean/SEM across participants...\n', thisExp);
    [grandMean_sound,  grandSEM_sound,  nSub_sound]  = grand_mean_sem(subjMeans_sound);
    [grandMean_type,   grandSEM_type,   nSub_type ]  = grand_mean_sem(subjMeans_type);
    [grandMean_effect, grandSEM_effect, nSub_eff  ]  = grand_mean_sem(subjMeans_effect);
    fprintf('[%s] Done.\n', thisExp);

    Nsub_title = subjCountLoaded;

    %% ---------- Save struct ----------
    GRAND.(thisExp).relTimes          = relTimes_ref;
    GRAND.(thisExp).soundOrder        = soundOrder;
    GRAND.(thisExp).typeOrder         = typeOrder;
    GRAND.(thisExp).effectOrder       = effectOrder;

    GRAND.(thisExp).subjMeans_sound   = subjMeans_sound;
    GRAND.(thisExp).subjMeans_type    = subjMeans_type;
    GRAND.(thisExp).subjMeans_effect  = subjMeans_effect;

    GRAND.(thisExp).grandMean_sound   = grandMean_sound;
    GRAND.(thisExp).grandSEM_sound    = grandSEM_sound;
    GRAND.(thisExp).nSub_sound        = nSub_sound;

    GRAND.(thisExp).grandMean_type    = grandMean_type;
    GRAND.(thisExp).grandSEM_type     = grandSEM_type;
    GRAND.(thisExp).nSub_type         = nSub_type;

    GRAND.(thisExp).grandMean_effect  = grandMean_effect;
    GRAND.(thisExp).grandSEM_effect   = grandSEM_effect;
    GRAND.(thisExp).nSub_effect       = nSub_eff;

    GRAND.(thisExp).Nsub_title        = Nsub_title;
    GRAND.(thisExp).Ntrials_sound     = totalTrials_sound;
    GRAND.(thisExp).Ntrials_effect    = totalTrials_effect;

    if saveGrandMat
        outMat = fullfile(grandOutFolder, sprintf('GRAND_%s_pupil.mat', thisExp));
        GRAND_this = GRAND.(thisExp); %#ok<NASGU>
        save(outMat, 'GRAND_this', '-v7.3');
        fprintf('[%s] Saved: %s\n', thisExp, outMat);
    end

    %% ===================== PLOTS =====================

    % --- Plot 1: LOCAL (Local Standard vs Local Deviant) ---
    fprintf('[%s] Plotting LOCAL (sound) grand average...\n', thisExp);
    fig1 = plot_grand_mean_sem_local( ...
        relTimes_ref, soundOrder, grandMean_sound, grandSEM_sound, ...
        thisExp, Nsub_title, totalTrials_sound);

    if saveFigures
        outPng = fullfile(grandOutFolder, sprintf('GRAND_%s_meanSEM_bySound.png', thisExp));
        exportgraphics(fig1, outPng, 'BackgroundColor', 'white', 'Resolution', 300);
        close(fig1);
        fprintf('[%s] Saved: %s\n', thisExp, outPng);
    end

    % --- Plot 2: TYPE (4 labels) ---
    fprintf('[%s] Plotting TYPE grand average...\n', thisExp);
    fig2 = plot_grand_mean_sem_by_type( ...
        relTimes_ref, typeOrder, grandMean_type, grandSEM_type, thisExp);

    if saveFigures
        outPng = fullfile(grandOutFolder, sprintf('GRAND_%s_meanSEM_byType.png', thisExp));
        exportgraphics(fig2, outPng, 'BackgroundColor', 'white', 'Resolution', 300);
        close(fig2);
        fprintf('[%s] Saved: %s\n', thisExp, outPng);
    end

    % --- Plot 3: ONE FIGURE, 2 PANELS (Local vs Global) ---
    fprintf('[%s] Plotting Local+Global (single figure)...\n', thisExp);
    fig3 = plot_local_global_onefig( ...
        relTimes_ref, ...
        grandMean_sound, grandSEM_sound, totalTrials_sound, ...
        grandMean_effect, grandSEM_effect, totalTrials_effect, ...
        thisExp, Nsub_title);

    if saveFigures
        outPng = fullfile(grandOutFolder, sprintf('GRAND_%s_meanSEM_LocalGlobal.png', thisExp));
        exportgraphics(fig3, outPng, 'BackgroundColor', 'white', 'Resolution', 300);
        close(fig3);
        fprintf('[%s] Saved: %s\n', thisExp, outPng);
    end
end

%% ===================== SAVE UPDATED T_LG =====================
save(infoMat,'T_LG','-append');
fprintf('\nT_LG updated and saved into: %s\n', infoMat);

excelFile = fullfile(outpath,'Info_EEG_LG.xlsx');
write_TLG_excel_exact(T_LG, excelFile);
fprintf('Excel export written: %s\n', excelFile);

fprintf('\nDONE. Outputs in: %s\n', grandOutFolder);


%% ===================== FUNCTIONS =====================

function T_LG = ensure_pupil_grandavg_columns(T_LG)
% ensure_pupil_grandavg_columns
% Adds (if missing) the per-subject pupil grand-average summary columns.
% IMPORTANT:
%   - Inserts them immediately AFTER 'BaselineAnchor_E2_Pupil'
%   - Orders columns as: ALL E1 first, then ALL E2
%   - Shifts any existing later columns to the right (nothing erased)
%
% Columns added (all numeric double):
%   E1:
%     Num_bad_trials_E1_Pupil
%     Num_good_trials_E1_Pupil
%     nStd_sound_E1_Pupil
%     nDev_sound_E1_Pupil
%     nLS_E1_Pupil   nLD_E1_Pupil   nGS_E1_Pupil   nGD_E1_Pupil   (Bekinschtein counts)
%
%   E2:
%     Num_bad_trials_E2_Pupil
%     Num_good_trials_E2_Pupil
%     nStd_sound_E2_Pupil
%     nDev_sound_E2_Pupil
%     nLS_E2_Pupil   nLD_E2_Pupil   nGS_E2_Pupil   nGD_E2_Pupil   (Bekinschtein counts)

anchor = 'BaselineAnchor_E2_Pupil';

if ~ismember(anchor, T_LG.Properties.VariableNames)
    warning('Column "%s" not found in T_LG. New columns will be appended at the end.', anchor);
end

% ---------- E1 block (ALL E1 first) ----------
newCols_E1 = { ...
    'Num_bad_trials_E1_Pupil', ...
    'Num_good_trials_E1_Pupil', ...
    'nStd_sound_E1_Pupil', ...
    'nDev_sound_E1_Pupil', ...
    'nLS_E1_Pupil', ...
    'nLD_E1_Pupil', ...
    'nGS_E1_Pupil', ...
    'nGD_E1_Pupil'};

% ---------- E2 block (then ALL E2) ----------
newCols_E2 = { ...
    'Num_bad_trials_E2_Pupil', ...
    'Num_good_trials_E2_Pupil', ...
    'nStd_sound_E2_Pupil', ...
    'nDev_sound_E2_Pupil', ...
    'nLS_E2_Pupil', ...
    'nLD_E2_Pupil', ...
    'nGS_E2_Pupil', ...
    'nGD_E2_Pupil'};

newCols = [newCols_E1, newCols_E2];

% ---------- Create missing columns (numeric) ----------
for k = 1:numel(newCols)
    vn = newCols{k};
    if ~ismember(vn, T_LG.Properties.VariableNames)
        T_LG.(vn) = nan(height(T_LG),1);
    end
end

% ---------- Move columns in ONE operation (preserves this order) ----------
if ismember(anchor, T_LG.Properties.VariableNames)
    T_LG = movevars(T_LG, newCols, 'After', anchor);
end

end



function m = subject_mean_curve(T, fieldName, condVar, condValue)
% subject-level mean curve (1 x nTime). If no trials -> NaNs.
    condVec = string(T.(condVar));
    idx = (condVec == string(condValue));

    nTime = numel(T.(fieldName){1});
    if ~any(idx)
        m = nan(1, nTime);
        return;
    end

    trialCells = T.(fieldName)(idx);
    pupMat = vertcat(trialCells{:}); % [nTrials x nTime]
    m = mean(pupMat, 1, 'omitnan');
end

function [grandMu, grandSem, nSub] = grand_mean_sem(subjMeansCell)
% equal-weight grand avg: mean over subject means, SEM over subjects
    nCond = numel(subjMeansCell);
    nTime = size(subjMeansCell{1}, 2);

    grandMu  = nan(nCond, nTime);
    grandSem = nan(nCond, nTime);
    nSub     = zeros(nCond, 1);

    for c = 1:nCond
        M = subjMeansCell{c};  % [nSub x nTime]

        contributed = any(~isnan(M), 2);
        nSub(c) = sum(contributed);
        if nSub(c) == 0
            continue;
        end

        grandMu(c, :) = mean(M, 1, 'omitnan');

        n_t  = sum(~isnan(M), 1);                    % per-timepoint N
        sd_t = std(M, 0, 1, 'omitnan');              % per-timepoint SD
        grandSem(c, :) = sd_t ./ sqrt(max(n_t, 1));  % per-timepoint SEM
    end
end

function fig = plot_grand_mean_sem_local(timeVector, soundOrder, grandMean, grandSEM, expLabel, Nsub_title, Ntrials_sound)
% Local Standard vs Local Deviant (clean legend: only line handles)

    % colors = [0 0.4470 0.7410;        % blue (Local Standard)
    %           0.8500 0.3250 0.0980];  % orange (Local Deviant)

    colors = [0 0.45 0.85;
    0.90 0.45 0.10];

    fig = figure('Units','normalized','OuterPosition',[0 0 1 1], ...
                 'Color','w', 'InvertHardcopy','off');
    ax = axes(fig); hold(ax,'on');
    set(ax,'Color','w','XColor','k','YColor','k');
    set(groot,'defaultAxesFontSize',11);  % affects all future axes
    axis tight

    title(ax, sprintf("Grand average | %s | Local | Nsub=%d | Ntrials LS=%d LD=%d", ...
        expLabel, Nsub_title, Ntrials_sound(1), Ntrials_sound(2)), FontSize=12);

    xlabel(ax,"Time (ms)", FontSize=12);
    ylabel(ax,"Pupil size (z-score)", FontSize=12);

    lineH = gobjects(numel(soundOrder),1);

    for s = 1:numel(soundOrder)
        mu  = grandMean(s, :);
        sem = grandSEM(s, :);

        if all(isnan(mu)), continue; end

        upper = mu + sem;
        lower = mu - sem;

        hFill = fill(ax, [timeVector fliplr(timeVector)], ...
                         [upper fliplr(lower)], ...
                         colors(s,:), 'FaceAlpha',0.15, 'EdgeColor','none');
        set(hFill,'HandleVisibility','off'); % IMPORTANT: do not pollute legend

        lineH(s) = plot(ax, timeVector, mu, ...
                        'LineWidth', 2, 'Color', colors(s,:), ...
                        'DisplayName', char(soundOrder(s)));
    end

    legend(ax, lineH(isgraphics(lineH)), 'Location','southwest', FontSize=12);
    grid(ax,'on');
    xline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');
    %hold(ax,'off');
end

function fig = plot_grand_mean_sem_by_type(timeVector, typeOrder, grandMean, grandSEM, expLabel)
% TYPE plot with correct legend mapping (only line handles)

    % colors = [ ...
    %     0.0000 0.4470 0.7410;  % Frequent Local Standard (blue)
    %     0.0000 0.4470 0.7410;  % Rare     Local Standard (blue)
    %     0.8500 0.3250 0.0980;  % Frequent Local Deviant  (orange)
    %     0.8500 0.3250 0.0980]; % Rare     Local Deviant  (orange)

    colors = [ ...
    0 0.45 0.85;
    0 0.45 0.85;
    0.90 0.45 0.10;
    0.90 0.45 0.10];

    lineStyles = {'--','-','--','-'};

    fig = figure('Units','normalized','OuterPosition',[0 0 1 1], ...
                 'Color','w', 'InvertHardcopy','off');
    ax = axes(fig); hold(ax,'on');
    set(ax,'Color','w','XColor','k','YColor','k');
    set(groot,'defaultAxesFontSize',11);  % affects all future axes
    axis tight

    title(ax, sprintf("Grand average | %s | Pupil respose by type", expLabel), FontSize=12);
    xlabel(ax,"Time (ms)", FontSize=12);
    ylabel(ax,"Pupil size (z-score)", FontSize=12);

    lineH = gobjects(numel(typeOrder),1);

    for k = 1:numel(typeOrder)
        mu  = grandMean(k, :);
        sem = grandSEM(k, :);

        if all(isnan(mu)), continue; end

        upper = mu + sem;
        lower = mu - sem;

        hFill = fill(ax, [timeVector fliplr(timeVector)], ...
                         [upper fliplr(lower)], ...
                         colors(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
        set(hFill,'HandleVisibility','off');

        lineH(k) = plot(ax, timeVector, mu, ...
                        'Color',     colors(k,:), ...
                        'LineStyle', lineStyles{k}, ...
                        'LineWidth', 2, ...
                        'DisplayName', char(typeOrder(k)));
    end

    legend(ax, lineH(isgraphics(lineH)), "Location","southwest", FontSize=12);
    grid(ax,'on');
    %hold(ax,'off');
    xline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');
end

function fig = plot_local_global_onefig(timeVector, mu_local, sem_local, Ntr_local, mu_global, sem_global, Ntr_global, expLabel, Nsub_title)
    % One figure, two stacked panels:
    %   Top    = Local Standard vs Local Deviant
    %   Bottom = Global Standard vs Global Deviant
    % Legend uses only the line handles. Titles contain Nsub and Ntrials per condition.

    % Colors: Standard=blue, Deviant=orange
    % cStd = [0 0.4470 0.7410];
    % cDev = [0.8500 0.3250 0.0980];
    cStd = [0 0.45 0.85];
    cDev = [0.90 0.45 0.10];

    fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.65], ...
             'InvertHardcopy','off');

    tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    % ---- TOP: Local ----
    ax1 = nexttile; hold(ax1,'on');
    set(ax1,'Color','w','XColor','k','YColor','k');
    set(groot,'defaultAxesFontSize',11);  % affects all future axes
    axis tight

    title(ax1, sprintf("Grand average | %s | Local | Nsub=%d | Ntrials LS=%d LD=%d", ...
        expLabel, Nsub_title, Ntr_local(1), Ntr_local(2)), FontSize=12);

    xlabel(ax1,"Time (ms)", FontSize=12);
    ylabel(ax1,"Pupil size (z-score)", FontSize=12);

    % LS (blue)
    upper = mu_local(1,:) + sem_local(1,:);
    lower = mu_local(1,:) - sem_local(1,:);
    hf = fill(ax1,[timeVector fliplr(timeVector)],[upper fliplr(lower)],cStd,'FaceAlpha',0.15,'EdgeColor','none');
    set(hf,'HandleVisibility','off');
    hLS = plot(ax1,timeVector,mu_local(1,:),'Color',cStd,'LineWidth',2,'DisplayName','Local Standard');

    % LD (orange)
    upper = mu_local(2,:) + sem_local(2,:);
    lower = mu_local(2,:) - sem_local(2,:);
    hf = fill(ax1,[timeVector fliplr(timeVector)],[upper fliplr(lower)],cDev,'FaceAlpha',0.15,'EdgeColor','none');
    set(hf,'HandleVisibility','off');
    hLD = plot(ax1,timeVector,mu_local(2,:),'Color',cDev,'LineWidth',2,'DisplayName','Local Deviant');

    xline(ax1,0,'k-','LineWidth',1,'HandleVisibility','off');
    legend(ax1,[hLS hLD],'Location','southwest', FontSize=12);
    grid(ax1,'on');
    %hold(ax1,'off');

    % ---- BOTTOM: Global ----
    ax2 = nexttile; hold(ax2,'on');
    set(ax2,'Color','w','XColor','k','YColor','k');

    title(ax2, sprintf("Grand average | %s | Global | Nsub=%d | Ntrials GS=%d GD=%d", ...
        expLabel, Nsub_title, Ntr_global(1), Ntr_global(2)), FontSize=12);

    xlabel(ax2,"Time (ms)", FontSize=12);
    ylabel(ax2,"Pupil size (z-score)", FontSize=12);

    % GS (blue)
    upper = mu_global(1,:) + sem_global(1,:);
    lower = mu_global(1,:) - sem_global(1,:);
    hf = fill(ax2,[timeVector fliplr(timeVector)],[upper fliplr(lower)],cStd,'FaceAlpha',0.15,'EdgeColor','none');
    set(hf,'HandleVisibility','off');
    hGS = plot(ax2,timeVector,mu_global(1,:),'Color',cStd,'LineWidth',2,'DisplayName','Global Standard');

    % GD (orange)
    upper = mu_global(2,:) + sem_global(2,:);
    lower = mu_global(2,:) - sem_global(2,:);
    hf = fill(ax2,[timeVector fliplr(timeVector)],[upper fliplr(lower)],cDev,'FaceAlpha',0.15,'EdgeColor','none');
    set(hf,'HandleVisibility','off');
    hGD = plot(ax2,timeVector,mu_global(2,:),'Color',cDev,'LineWidth',2,'DisplayName','Global Deviant');

    xline(ax2,0,'k-','LineWidth',1,'HandleVisibility','off');
    legend(ax2,[hGS hGD],'Location','southwest', FontSize=12);
    grid(ax2,'on');
    %hold(ax2,'off');
end
