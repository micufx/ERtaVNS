clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% D3_Pupil_ERP_ROI_correlation.m
%
% PURPOSE
%   Participant-level correlations between pupil deviance effects and
%   centro-parietal ROI ERP deviance effects in the local-global paradigm.
%
%   This is an ROI adaptation of C5_Pupil_ERP_correlation.m. The only
%   substantive change is the EEG feature:
%
%       OLD: one predefined channel (Pz)
%       NEW: E1-derived centro-parietal ROI, separately for Local and Global
%
% E1-DERIVED ROI DEFINITION
%   The script reconstructs the E1 subject-level Local and Global ERPs and
%   reproduces the group logic used in B7_EEG_stats.m:
%
%       1. Paired Dev - Std t-map across subjects, 0-700 ms
%       2. max-|t| sign-flip permutation FWE correction
%       3. significant time requires:
%              corrected p <= 0.01
%              at least 5 electrodes
%              at least 20 consecutive samples
%       4. Within the surviving E1 time samples, channels are retained when
%          they show corrected p <= 0.01 and belong to the predefined
%          centro-parietal channel pool.
%       5. The Local and Global ROIs are therefore allowed to contain
%          different channels and different E1-derived time windows.
%
%   The selected E1 ROI and time window are then applied unchanged to BOTH
%   E1 and E2. This follows the user's requested logic and makes E2 the more
%   independent transfer test.
%
% IMPORTANT CIRCULARITY NOTE
%   For E1, the same participants define the significant Dev-Std ROI/window
%   and are then used in the E1 pupil-ERP correlation. This is exploratory
%   and is more circular than Alamia et al., whose initial EEG clusters were
%   defined from responses versus baseline rather than from the rare-
%   frequent contrast itself. The E2 application is less circular because
%   the E1-derived ROI/window is transferred unchanged to E2.
%
% FEATURES
%   For each participant and experiment:
%       LOCAL pupil effect  = Local Deviant  - Local Standard
%       GLOBAL pupil effect = Global Deviant - Global Standard
%       LOCAL ROI ERP effect  = ROI mean(Local Deviant - Local Standard)
%       GLOBAL ROI ERP effect = ROI mean(Global Deviant - Global Standard)
%
% PRIMARY / ALAMIA-CLOSEST ANALYSIS
%       Pupil signed absolute peak effect
%           versus
%       ERP mean effect averaged across the E1-derived ROI and time window
%
% ADDITIONAL ANALYSES
%       Pupil peak vs ERP peak
%       Pupil mean vs ERP mean
%       Local, Global, and exploratory General composite
%
% OUTPUT
%   - E1_ROI_definition_channels.csv
%   - E1_ROI_definition_windows.csv
%   - E1_ROI_definition.mat
%   - Pupil_ERP_ROI_subject_features.csv
%   - Pupil_ERP_ROI_correlation_results.csv
%   - Pupil_ERP_ROI_correlation_results.mat
%   - Six figures matching the previous script, now using ROI ERP features
%
% Miguel Contreras-Altamirano, 2026

%% ===================== USER SETTINGS =====================

cfg = struct();

% Main data paths
cfg.eeglabPath = 'L:\Downloads\eeglab2023.0\';
cfg.outroot    = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';
cfg.outFolder  = fullfile(cfg.outroot, 'GRAND_AVG', 'Pupil_ERP_ROI_correlation');

if ~exist(cfg.outFolder, 'dir')
    mkdir(cfg.outFolder);
end

% Data definitions
cfg.pupilField = 'PUPIL_Z';

cfg.targets = [ ...
    "Frequent Local Standard", ...
    "Rare Local Standard", ...
    "Frequent Local Deviant", ...
    "Rare Local Deviant"];

% Exclusions used in the existing pupil and EEG scripts
cfg.excludeSubjects = ["S21","S12","S23","S02","S29"];

% Pupil significant group windows from the manuscript
cfg.pupilWindows.E1.Local  = [44 1232];
cfg.pupilWindows.E1.Global = [288 2000];
cfg.pupilWindows.E2.Local  = [0 690];
cfg.pupilWindows.E2.Global = [446 1818];

% -------- E1 EEG ROI derivation: matches B7_EEG_stats.m --------
cfg.roi.statsWindowMS = [0 700];
cfg.roi.alphaCorrected = 0.01;
cfg.roi.minElectrodes   = 5;
cfg.roi.minConsecutiveSamples = 20;

cfg.roi.fweNPerm       = 5000;
cfg.roi.fweSeed        = 777;
cfg.roi.verbose        = true;
cfg.roi.progressStep   = 200;

% A centro-parietal anatomical pool. Only E1-significant red-dot channels
% that also belong to this list are retained. Edit this list only if you
% want a narrower or wider anatomical definition.
cfg.roi.centroParietalPool = [ ...
    "C3","Cz","C4", ...
    "CP5","CP1","CP2","CP6", ...
    "P3","Pz","P4","P7","P8"];

% A candidate channel must be significant for at least this many samples
% inside the surviving E1 extent-mask interval.
cfg.roi.minSignificantSamplesPerChannel = 1;

% Peak definition for difference waves:
%   'signed_abs_max' = sample with largest absolute magnitude, retaining sign
%   'positive_max'   = largest positive sample in the window
cfg.peakDefinition = 'signed_abs_max';

% Correlation
cfg.corrTypePrimary = 'Pearson';
cfg.nPerm           = 10000;
cfg.permutationSeed = 2026;

% Bayesian analysis
cfg.doBayesianCorrelation = true;

% Optional explicit Rscript path, for example:
% cfg.rscriptPath = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe';
% Leave empty to use "Rscript" from the system PATH.
cfg.rscriptPath = '';

cfg.bayesPosteriorIterations = 30000;
cfg.requireExactBayesFactor  = false;

% Figures
cfg.makeAlamiaLikeFigures = true;
cfg.saveFigures           = true;
cfg.showSubjectLabels     = false;
cfg.exportDPI             = 300;

% Requested colors
cfg.colors = [ ...
    0.00 0.45 0.85; ...  % participant dots: blue
    0.90 0.45 0.10];     % regression line: orange

% General composite
% 'raw_equal_weight' or 'zscore_equal_weight'
cfg.generalCompositeMode = 'raw_equal_weight';

%% ===================== START EEGLAB =====================

if exist('pop_loadset', 'file') ~= 2
    if exist(cfg.eeglabPath, 'dir')
        addpath(cfg.eeglabPath);
    else
        error('EEGLAB path not found:\n%s', cfg.eeglabPath);
    end
end

if exist('pop_loadset', 'file') ~= 2
    try
        eeglab('nogui');
    catch
        error('Could not initialize EEGLAB.');
    end
end

%% ===================== FIND PARTICIPANTS =====================

subjects = list_subjects(cfg.outroot, cfg.excludeSubjects);

if isempty(subjects)
    error('No participant folders were found under:\n%s', cfg.outroot);
end

fprintf('\nParticipants found after exclusions: %d\n', numel(subjects));
fprintf('%s\n', strjoin(cellstr(subjects), ', '));

%% ===================== DERIVE E1 CENTRO-PARIETAL ROIs =====================

fprintf('\n============================================================\n');
fprintf('DERIVING E1 SIGNIFICANT CENTRO-PARIETAL ROIs\n');
fprintf('============================================================\n');

[ROI, ROI_DIAGNOSTICS] = derive_e1_significant_rois(subjects, cfg);

% The E1-derived windows are reused unchanged for E2.
cfg.erpWindows.Local  = ROI.Local.WindowMS;
cfg.erpWindows.Global = ROI.Global.WindowMS;

fprintf('\nSelected LOCAL ROI channels: %s\n', ...
    strjoin(cellstr(ROI.Local.ChannelLabels), ', '));
fprintf('Selected LOCAL E1 window: %.0f to %.0f ms\n', ...
    ROI.Local.WindowMS(1), ROI.Local.WindowMS(2));

fprintf('\nSelected GLOBAL ROI channels: %s\n', ...
    strjoin(cellstr(ROI.Global.ChannelLabels), ', '));
fprintf('Selected GLOBAL E1 window: %.0f to %.0f ms\n', ...
    ROI.Global.WindowMS(1), ROI.Global.WindowMS(2));

fprintf(['\nCAUTION: E1 ROI/window selection and the E1 correlation use the same ' ...
    'participants. Treat E1 as exploratory. E2 uses the transferred E1 ROI/window.\n']);

save_roi_definition_outputs(ROI, ROI_DIAGNOSTICS, cfg);

%% ===================== EXTRACT SUBJECT FEATURES =====================

FEATURES = table();
expList = ["E1","E2"];

for e = 1:numel(expList)

    expKey = expList(e);

    fprintf('\n============================================================\n');
    fprintf('EXTRACTING PUPIL AND ROI ERP FEATURES: %s\n', expKey);
    fprintf('============================================================\n');

    for s = 1:numel(subjects)

        subject = subjects(s);

        pupilFile = fullfile(cfg.outroot, char(subject), ...
            sprintf('%s_%s_pupil_preproc.mat', subject, expKey));

        eegFile = fullfile(cfg.outroot, char(subject), ...
            sprintf('auditory_all_LG_%s_%s.set', subject, expKey));

        if exist(pupilFile, 'file') ~= 2
            fprintf('  [%s %s] Missing pupil file -> skipped\n', subject, expKey);
            continue;
        end

        if exist(eegFile, 'file') ~= 2
            fprintf('  [%s %s] Missing EEG file -> skipped\n', subject, expKey);
            continue;
        end

        try
            P = extract_pupil_subject_features(pupilFile, expKey, cfg);
        catch
            warning('[%s %s] Pupil extraction failed.', subject, expKey);
            continue;
        end

        try
            E = extract_erp_subject_features(eegFile, expKey, cfg, ROI);
        catch
            warning('[%s %s] ERP extraction failed.', subject, expKey);
            continue;
        end

        row = struct();

        row.Subject    = string(subject);
        row.Experiment = string(expKey);

        % Pupil local
        row.PupilLocalStandardMean = P.Local.StandardMean;
        row.PupilLocalDeviantMean  = P.Local.DeviantMean;
        row.PupilMeanLocal         = P.Local.MeanEffect;
        row.PupilPeakLocal         = P.Local.PeakEffect;
        row.PupilPeakLatencyLocal  = P.Local.PeakLatency;

        % Pupil global
        row.PupilGlobalStandardMean = P.Global.StandardMean;
        row.PupilGlobalDeviantMean  = P.Global.DeviantMean;
        row.PupilMeanGlobal         = P.Global.MeanEffect;
        row.PupilPeakGlobal         = P.Global.PeakEffect;
        row.PupilPeakLatencyGlobal  = P.Global.PeakLatency;

        % ERP local
        row.ERPLocalStandardMean = E.Local.StandardMean;
        row.ERPLocalDeviantMean  = E.Local.DeviantMean;
        row.ERPMeanLocal         = E.Local.MeanEffect;
        row.ERPPeakLocal         = E.Local.PeakEffect;
        row.ERPPeakLatencyLocal  = E.Local.PeakLatency;

        % ERP global
        row.ERPGlobalStandardMean = E.Global.StandardMean;
        row.ERPGlobalDeviantMean  = E.Global.DeviantMean;
        row.ERPMeanGlobal         = E.Global.MeanEffect;
        row.ERPPeakGlobal         = E.Global.PeakEffect;
        row.ERPPeakLatencyGlobal  = E.Global.PeakLatency;

        row.PupilLocalWindowStart  = cfg.pupilWindows.(char(expKey)).Local(1);
        row.PupilLocalWindowEnd    = cfg.pupilWindows.(char(expKey)).Local(2);
        row.PupilGlobalWindowStart = cfg.pupilWindows.(char(expKey)).Global(1);
        row.PupilGlobalWindowEnd   = cfg.pupilWindows.(char(expKey)).Global(2);

        row.ERPLocalWindowStart  = cfg.erpWindows.Local(1);
        row.ERPLocalWindowEnd    = cfg.erpWindows.Local(2);
        row.ERPGlobalWindowStart = cfg.erpWindows.Global(1);
        row.ERPGlobalWindowEnd   = cfg.erpWindows.Global(2);

        row.ERPLocalROINChannels  = numel(ROI.Local.ChannelLabels);
        row.ERPGlobalROINChannels = numel(ROI.Global.ChannelLabels);
        row.ERPLocalROIChannels   = strjoin(ROI.Local.ChannelLabels, ';');
        row.ERPGlobalROIChannels  = strjoin(ROI.Global.ChannelLabels, ';');

        FEATURES = [FEATURES; struct2table(row, 'AsArray', true)]; %#ok<AGROW>

        fprintf('  [%s %s] completed\n', subject, expKey);

    end

end

if isempty(FEATURES)
    error('No complete pupil + EEG participant records were extracted.');
end

%% ===================== COMPUTE GENERAL COMPOSITES =====================

FEATURES = add_general_composites(FEATURES, cfg.generalCompositeMode);

% Sort by experiment and subject number
FEATURES.SubjectNumber = subject_number(FEATURES.Subject);
FEATURES = sortrows(FEATURES, {'Experiment','SubjectNumber'});

%% ===================== SAVE EXTRACTED FEATURES =====================

featureCsv = fullfile(cfg.outFolder, 'Pupil_ERP_ROI_subject_features.csv');
writetable(FEATURES, featureCsv);

fprintf('\nSaved subject-level features:\n  %s\n', featureCsv);

%% ===================== DEFINE CORRELATIONS =====================

% Each correlation is computed separately within E1 and E2.
%
% "Alamia" measure:
%   Pupil peak effect versus ERP mean effect.
%
% Requested matched measures:
%   Peak: pupil peak effect versus ERP peak effect.
%   Mean: pupil mean effect versus ERP mean effect.

ANALYSES = build_analysis_definitions();

%% ===================== RUN CORRELATIONS =====================

RESULTS = table();
BAYES_LONG = table();

for a = 1:height(ANALYSES)

    expKey = ANALYSES.Experiment(a);
    xVar   = ANALYSES.PupilVariable(a);
    yVar   = ANALYSES.ERPVariable(a);

    idx = FEATURES.Experiment == expKey;

    x = FEATURES.(char(xVar))(idx);
    y = FEATURES.(char(yVar))(idx);
    subj = FEATURES.Subject(idx);

    ST = run_correlation(x, y, subj, cfg, a);

    resultRow = table( ...
        ANALYSES.AnalysisID(a), ...
        ANALYSES.AnalysisFamily(a), ...
        expKey, ...
        ANALYSES.Measure(a), ...
        ANALYSES.Contrast(a), ...
        xVar, ...
        yVar, ...
        ST.N, ...
        ST.rPearson, ...
        ST.pPearson, ...
        ST.pPermutation, ...
        ST.rSpearman, ...
        ST.pSpearman, ...
        ST.CI95_Pearson(1), ...
        ST.CI95_Pearson(2), ...
        'VariableNames', { ...
        'AnalysisID','AnalysisFamily','Experiment','Measure','Contrast', ...
        'PupilVariable','ERPVariable','N','rPearson','pPearson', ...
        'pPermutation','rhoSpearman','pSpearman', ...
        'PearsonCI95_Low','PearsonCI95_High'});

    RESULTS = [RESULTS; resultRow]; %#ok<AGROW>

    if ST.N >= 4
        tmp = table( ...
            repmat(ANALYSES.AnalysisID(a), ST.N, 1), ...
            ST.subjects, ST.x, ST.y, ...
            'VariableNames', {'AnalysisID','Subject','Pupil','ERP'});

        BAYES_LONG = [BAYES_LONG; tmp]; %#ok<AGROW>
    end

end

%% ===================== BAYESIAN CORRELATIONS =====================

RESULTS.BF10 = nan(height(RESULTS),1);
RESULTS.BFmethod = strings(height(RESULTS),1);
RESULTS.BayesCrI95_Low  = nan(height(RESULTS),1);
RESULTS.BayesCrI95_High = nan(height(RESULTS),1);

if cfg.doBayesianCorrelation

    [BAYES_RESULTS, exactBayesWorked, bayesMessage] = ...
        run_bayesian_correlations_batch(BAYES_LONG, cfg);

    fprintf('\nBayesian analysis: %s\n', bayesMessage);

    if exactBayesWorked

        [tf, loc] = ismember(RESULTS.AnalysisID, BAYES_RESULTS.AnalysisID);

        RESULTS.BF10(tf) = BAYES_RESULTS.BF10(loc(tf));
        RESULTS.BFmethod(tf) = BAYES_RESULTS.Method(loc(tf));
        RESULTS.BayesCrI95_Low(tf)  = BAYES_RESULTS.CrI95_Low(loc(tf));
        RESULTS.BayesCrI95_High(tf) = BAYES_RESULTS.CrI95_High(loc(tf));

    else

        if cfg.requireExactBayesFactor
            error(['Exact R-based Bayesian correlations could not be computed.\n' ...
                   'Reason: %s'], bayesMessage);
        end

        for a = 1:height(RESULTS)
            RESULTS.BF10(a) = bic_bf10_from_correlation( ...
                RESULTS.rPearson(a), RESULTS.N(a));
            RESULTS.BFmethod(a) = "BIC approximation";
        end

    end

else

    RESULTS.BFmethod(:) = "Not computed";

end

%% ===================== MULTIPLE-COMPARISON CORRECTION =====================

% Holm correction is applied separately within each analysis family:
%   Alamia-like, Peak-Peak, and Mean-Mean.
%
% This preserves the distinction between primary/paper-close and
% exploratory analyses.

RESULTS.pPearsonHolm    = nan(height(RESULTS),1);
RESULTS.pPermutationHolm = nan(height(RESULTS),1);

families = unique(RESULTS.AnalysisFamily);

for f = 1:numel(families)
    idx = RESULTS.AnalysisFamily == families(f);
    RESULTS.pPearsonHolm(idx) = holm_adjust(RESULTS.pPearson(idx));
    RESULTS.pPermutationHolm(idx) = holm_adjust(RESULTS.pPermutation(idx));
end

%% ===================== SAVE RESULTS =====================

resultsCsv = fullfile(cfg.outFolder, 'Pupil_ERP_ROI_correlation_results.csv');
writetable(RESULTS, resultsCsv);

resultsMat = fullfile(cfg.outFolder, 'Pupil_ERP_ROI_correlation_results.mat');
save(resultsMat, 'FEATURES', 'RESULTS', 'ANALYSES', 'ROI', 'ROI_DIAGNOSTICS', 'cfg', '-v7.3');

fprintf('Saved correlation table:\n  %s\n', resultsCsv);
fprintf('Saved MATLAB results:\n  %s\n', resultsMat);

%% ===================== PRINT COMPACT SUMMARY =====================

fprintf('\n============================================================\n');
fprintf('PUPIL-ERP CORRELATION SUMMARY\n');
fprintf('============================================================\n');

for a = 1:height(RESULTS)

    fprintf(['%s | %s | %s | N=%d | r=%.3f | p=%.4f | ' ...
             'perm p=%.4f | BF10=%s\n'], ...
        RESULTS.AnalysisFamily(a), ...
        RESULTS.Experiment(a), ...
        RESULTS.Contrast(a), ...
        RESULTS.N(a), ...
        RESULTS.rPearson(a), ...
        RESULTS.pPearson(a), ...
        RESULTS.pPermutation(a), ...
        format_bf(RESULTS.BF10(a)));

end

%% ===================== CREATE FIGURES =====================

if cfg.makeAlamiaLikeFigures

    % 1. Paper-closest: pupil peak versus ERP mean, local/global
    fig1 = plot_local_global_figure( ...
        FEATURES, RESULTS, ...
        "Alamia", ...
        'Pupil peak effect vs centro-parietal ROI ERP mean effect', ...
        cfg);

    save_figure_pair(fig1, cfg.outFolder, ...
        'Fig_ROI_PupilPeak_vs_ERPMean_LOCAL_GLOBAL_AlamiaLike', cfg);

    % 2. Paper-closest: pupil peak versus ERP mean, general composite
    fig2 = plot_general_figure( ...
        FEATURES, RESULTS, ...
        "Alamia", ...
        'Pupil peak effect vs centro-parietal ROI ERP mean effect: general composite', ...
        cfg);

    save_figure_pair(fig2, cfg.outFolder, ...
        'Fig_ROI_PupilPeak_vs_ERPMean_GENERAL_AlamiaLike', cfg);

end

% 3. Requested peak-peak local/global
fig3 = plot_local_global_figure( ...
    FEATURES, RESULTS, ...
    "PeakPeak", ...
    'Pupil peak effect vs centro-parietal ROI ERP peak effect', ...
    cfg);

save_figure_pair(fig3, cfg.outFolder, ...
    'Fig_ROI_PupilPeak_vs_ERPPeak_LOCAL_GLOBAL', cfg);

% 4. Requested peak-peak general
fig4 = plot_general_figure( ...
    FEATURES, RESULTS, ...
    "PeakPeak", ...
    'Pupil peak effect vs centro-parietal ROI ERP peak effect: general composite', ...
    cfg);

save_figure_pair(fig4, cfg.outFolder, ...
    'Fig_ROI_PupilPeak_vs_ERPPeak_GENERAL', cfg);

% 5. Requested mean-mean local/global
fig5 = plot_local_global_figure( ...
    FEATURES, RESULTS, ...
    "MeanMean", ...
    'Pupil mean effect vs centro-parietal ROI ERP mean effect', ...
    cfg);

save_figure_pair(fig5, cfg.outFolder, ...
    'Fig_ROI_PupilMean_vs_ERPMean_LOCAL_GLOBAL', cfg);

% 6. Requested mean-mean general
fig6 = plot_general_figure( ...
    FEATURES, RESULTS, ...
    "MeanMean", ...
    'Pupil mean effect vs centro-parietal ROI ERP mean effect: general composite', ...
    cfg);

save_figure_pair(fig6, cfg.outFolder, ...
    'Fig_ROI_PupilMean_vs_ERPMean_GENERAL', cfg);

fprintf('\nDone. All outputs were saved in:\n  %s\n', cfg.outFolder);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function subjects = list_subjects(outroot, excludeSubjects)

    d = dir(fullfile(outroot, 'S*'));
    d = d([d.isdir]);

    subjects = string({d.name});
    subjects = subjects(~ismember(subjects, [".","..","GRAND_AVG"]));

    keep = ~cellfun('isempty', regexp(cellstr(subjects), '^S\d+$', 'once'));
    subjects = subjects(keep);

    if ~isempty(excludeSubjects)
        subjects = subjects(~ismember(lower(subjects), lower(string(excludeSubjects))));
    end

    nums = subject_number(subjects);
    [~, order] = sort(nums);
    subjects = subjects(order);

end

function nums = subject_number(subjects)

    subjects = string(subjects);
    nums = nan(numel(subjects),1);

    for i = 1:numel(subjects)
        tok = regexp(char(subjects(i)), '\d+', 'match', 'once');
        if ~isempty(tok)
            nums(i) = str2double(tok);
        end
    end

end

function P = extract_pupil_subject_features(pupilFile, expKey, cfg)

    S = load(pupilFile, 'dataTableClean', 'relTimes');

    if ~isfield(S, 'dataTableClean') || isempty(S.dataTableClean)
        error('dataTableClean is missing or empty.');
    end

    if ~isfield(S, 'relTimes') || isempty(S.relTimes)
        error('relTimes is missing or empty.');
    end

    T = S.dataTableClean;
    times = double(S.relTimes(:))';

    if ~ismember(cfg.pupilField, T.Properties.VariableNames)
        error('Pupil field "%s" is missing.', cfg.pupilField);
    end

    if ~ismember('sound', T.Properties.VariableNames)
        error('Column "sound" is missing from dataTableClean.');
    end

    if ~ismember('type', T.Properties.VariableNames)
        error('Column "type" is missing from dataTableClean.');
    end

    soundVec = string(T.sound);
    typeVec  = string(T.type);

    % Local pupil conditions follow the current pupil script directly:
    % all Local Standard trials versus all Local Deviant trials.
    mLS = pupil_mean_curve(T, cfg.pupilField, soundVec == "Local Standard");
    mLD = pupil_mean_curve(T, cfg.pupilField, soundVec == "Local Deviant");

    % Global pupil conditions:
    % Frequent = Global Standard; Rare = Global Deviant.
    isGS = startsWith(typeVec, "Frequent", 'IgnoreCase', true);
    isGD = startsWith(typeVec, "Rare", 'IgnoreCase', true);

    mGS = pupil_mean_curve(T, cfg.pupilField, isGS);
    mGD = pupil_mean_curve(T, cfg.pupilField, isGD);

    if numel(mLS) ~= numel(times) || numel(mLD) ~= numel(times) || ...
            numel(mGS) ~= numel(times) || numel(mGD) ~= numel(times)
        error('Pupil curve length does not match relTimes.');
    end

    localWindow  = cfg.pupilWindows.(char(expKey)).Local;
    globalWindow = cfg.pupilWindows.(char(expKey)).Global;

    P.Local = effect_features_from_curves( ...
        mLS, mLD, times, localWindow, cfg.peakDefinition);

    P.Global = effect_features_from_curves( ...
        mGS, mGD, times, globalWindow, cfg.peakDefinition);

end

function curve = pupil_mean_curve(T, fieldName, trialMask)

    trialMask = logical(trialMask(:));

    if numel(trialMask) ~= height(T)
        error('Pupil trial mask length does not match table height.');
    end

    if ~any(trialMask)
        error('No trials were found for one pupil condition.');
    end

    cells = T.(fieldName)(trialMask);

    if iscell(cells)
        mat = vertcat(cells{:});
    elseif isnumeric(cells)
        mat = double(cells);
    else
        error('Unsupported pupil data type in %s.', fieldName);
    end

    curve = mean(double(mat), 1, 'omitnan');

end

function [ROI, D] = derive_e1_significant_rois(subjects, cfg)
% Derive Local and Global E1 centro-parietal ROIs from corrected group maps.

    E1_Lstd = [];
    E1_Ldev = [];
    E1_Gstd = [];
    E1_Gdev = [];
    usedSubjects = strings(0,1);

    templateLabels = strings(0,1);
    templateTimes = [];

    for s = 1:numel(subjects)

        subject = subjects(s);
        eegFile = fullfile(cfg.outroot, char(subject), ...
            sprintf('auditory_all_LG_%s_E1.set', subject));

        if exist(eegFile, 'file') ~= 2
            fprintf('  [ROI E1 %s] EEG file missing -> skipped\n', subject);
            continue;
        end

        [filePath, fileName, fileExt] = fileparts(eegFile);
        EEG = pop_loadset('filename', [fileName fileExt], 'filepath', filePath);
        EEG = eeg_checkset(EEG);

        if isempty(templateLabels)
            templateLabels = string({EEG.chanlocs.labels});
            templateTimes = double(EEG.times(:))';
        else
            EEG = reorder_eeg_to_template(EEG, templateLabels, templateTimes);
        end

        try
            [Lstd, Ldev, Gstd, Gdev] = compute_subject_lg_erps(EEG, cfg.targets);
        catch
            warning('[ROI E1 %s] Could not reconstruct conditions.', subject);
            continue;
        end

        E1_Lstd(end+1,:,:) = Lstd; %#ok<AGROW>
        E1_Ldev(end+1,:,:) = Ldev; %#ok<AGROW>
        E1_Gstd(end+1,:,:) = Gstd; %#ok<AGROW>
        E1_Gdev(end+1,:,:) = Gdev; %#ok<AGROW>
        usedSubjects(end+1,1) = subject; %#ok<AGROW>

        fprintf('  [ROI E1 %s] included\n', subject);

    end

    nSub = size(E1_Lstd,1);
    if nSub < 4
        error('Only %d E1 participants were available for ROI derivation.', nSub);
    end

    statIdx = templateTimes >= cfg.roi.statsWindowMS(1) & ...
              templateTimes <= cfg.roi.statsWindowMS(2);

    statTimes = templateTimes(statIdx);

    LstdStat = E1_Lstd(:,:,statIdx);
    LdevStat = E1_Ldev(:,:,statIdx);
    GstdStat = E1_Gstd(:,:,statIdx);
    GdevStat = E1_Gdev(:,:,statIdx);

    fprintf('\nComputing E1 Local max-|t| FWE p-map (%d permutations)...\n', ...
        cfg.roi.fweNPerm);
    [pLocFWE, infoLoc] = fwe_maxT_paired_roi( ...
        LdevStat, LstdStat, cfg.roi.fweNPerm, cfg.roi.fweSeed, ...
        cfg.roi.verbose, cfg.roi.progressStep);

    fprintf('\nComputing E1 Global max-|t| FWE p-map (%d permutations)...\n', ...
        cfg.roi.fweNPerm);
    [pGloFWE, infoGlo] = fwe_maxT_paired_roi( ...
        GdevStat, GstdStat, cfg.roi.fweNPerm, cfg.roi.fweSeed, ...
        cfg.roi.verbose, cfg.roi.progressStep);

    maskLoc = triple_criterion_time_mask_roi( ...
        pLocFWE, cfg.roi.alphaCorrected, ...
        cfg.roi.minElectrodes, cfg.roi.minConsecutiveSamples);

    maskGlo = triple_criterion_time_mask_roi( ...
        pGloFWE, cfg.roi.alphaCorrected, ...
        cfg.roi.minElectrodes, cfg.roi.minConsecutiveSamples);

    ROI = struct();

    ROI.Local = select_effect_roi( ...
        'Local', pLocFWE, maskLoc, statTimes, templateLabels, cfg);

    ROI.Global = select_effect_roi( ...
        'Global', pGloFWE, maskGlo, statTimes, templateLabels, cfg);

    ROI.SelectionDataset = 'E1';
    ROI.SelectionN = nSub;
    ROI.SelectionSubjects = usedSubjects;
    ROI.CentroParietalPool = cfg.roi.centroParietalPool;
    ROI.AlphaCorrected = cfg.roi.alphaCorrected;
    ROI.MinElectrodes = cfg.roi.minElectrodes;
    ROI.MinConsecutiveSamples = cfg.roi.minConsecutiveSamples;
    ROI.FWENPerm = cfg.roi.fweNPerm;
    ROI.FWESeed = cfg.roi.fweSeed;

    D = struct();
    D.Subjects = usedSubjects;
    D.N = nSub;
    D.ChannelLabels = templateLabels;
    D.Times = statTimes;
    D.Local.pFWE = pLocFWE;
    D.Local.TimeMask = maskLoc;
    D.Local.Info = infoLoc;
    D.Local.Lstd = LstdStat;
    D.Local.Ldev = LdevStat;
    D.Global.pFWE = pGloFWE;
    D.Global.TimeMask = maskGlo;
    D.Global.Info = infoGlo;
    D.Global.Gstd = GstdStat;
    D.Global.Gdev = GdevStat;

end

function R = select_effect_roi(effectName, pFWE, timeMask, times, labels, cfg)

    if ~any(timeMask)
        error(['No E1 %s time interval survived the corrected p-map plus ' ...
            'extent criterion. The requested ROI cannot be defined.'], effectName);
    end

    segments = logical_runs_to_windows(timeMask, times);

    % User requested the significant overall time window.
    overallWindow = [times(find(timeMask,1,'first')), ...
                     times(find(timeMask,1,'last'))];

    if size(segments,1) > 1
        warning(['E1 %s contained %d separated significant segments. ' ...
            'The ROI feature uses the requested overall first-to-last window ' ...
            '[%.0f %.0f] ms. Segment boundaries are saved separately.'], ...
            effectName, size(segments,1), overallWindow(1), overallWindow(2));
    end

    candidateIdx = labels_to_indices(labels, cfg.roi.centroParietalPool);

    if isempty(candidateIdx)
        error('None of the configured centro-parietal channels exist in the EEG montage.');
    end

    sigMapWithinExtent = (pFWE <= cfg.roi.alphaCorrected) & timeMask;
    sigCounts = sum(sigMapWithinExtent, 2);

    selectedIdx = candidateIdx( ...
        sigCounts(candidateIdx) >= cfg.roi.minSignificantSamplesPerChannel);

    if isempty(selectedIdx)
        candidateText = strjoin(cellstr(labels(candidateIdx)), ', ');
        error(['No centro-parietal channel was significantly marked within ' ...
            'the E1 %s extent interval. Candidate pool present in data: %s'], ...
            effectName, candidateText);
    end

    R = struct();
    R.Effect = string(effectName);
    R.ChannelLabels = labels(selectedIdx);
    R.ChannelIndicesInTemplate = selectedIdx(:);
    R.SignificantSampleCounts = sigCounts(selectedIdx);
    R.WindowMS = overallWindow;
    R.SignificantSegmentsMS = segments;
    R.TimeMask = timeMask;
    R.SelectedChannelMask = false(numel(labels),1);
    R.SelectedChannelMask(selectedIdx) = true;

end

function windows = logical_runs_to_windows(mask, times)

    mask = logical(mask(:)');
    times = double(times(:)');

    d = diff([false mask false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    windows = nan(numel(starts),2);
    for k = 1:numel(starts)
        windows(k,:) = [times(starts(k)), times(ends(k))];
    end

end

function save_roi_definition_outputs(ROI, D, cfg)

    effect = strings(0,1);
    channel = strings(0,1);
    significantSamples = [];
    windowStart = [];
    windowEnd = [];

    effects = ["Local","Global"];

    for e = 1:numel(effects)
        key = char(effects(e));
        labs = ROI.(key).ChannelLabels;
        counts = ROI.(key).SignificantSampleCounts;

        for c = 1:numel(labs)
            effect(end+1,1) = effects(e); %#ok<AGROW>
            channel(end+1,1) = labs(c); %#ok<AGROW>
            significantSamples(end+1,1) = counts(c); %#ok<AGROW>
            windowStart(end+1,1) = ROI.(key).WindowMS(1); %#ok<AGROW>
            windowEnd(end+1,1) = ROI.(key).WindowMS(2); %#ok<AGROW>
        end
    end

    Tchannels = table(effect, channel, significantSamples, windowStart, windowEnd, ...
        'VariableNames', {'Effect','Channel','SignificantSamples', ...
        'OverallWindowStartMS','OverallWindowEndMS'});

    writetable(Tchannels, fullfile(cfg.outFolder, 'E1_ROI_definition_channels.csv'));

    effectW = strings(0,1);
    segment = [];
    segmentStart = [];
    segmentEnd = [];
    overallStart = [];
    overallEnd = [];

    for e = 1:numel(effects)
        key = char(effects(e));
        seg = ROI.(key).SignificantSegmentsMS;

        for k = 1:size(seg,1)
            effectW(end+1,1) = effects(e); %#ok<AGROW>
            segment(end+1,1) = k; %#ok<AGROW>
            segmentStart(end+1,1) = seg(k,1); %#ok<AGROW>
            segmentEnd(end+1,1) = seg(k,2); %#ok<AGROW>
            overallStart(end+1,1) = ROI.(key).WindowMS(1); %#ok<AGROW>
            overallEnd(end+1,1) = ROI.(key).WindowMS(2); %#ok<AGROW>
        end
    end

    Twindows = table(effectW, segment, segmentStart, segmentEnd, overallStart, overallEnd, ...
        'VariableNames', {'Effect','Segment','SegmentStartMS','SegmentEndMS', ...
        'OverallWindowStartMS','OverallWindowEndMS'});

    writetable(Twindows, fullfile(cfg.outFolder, 'E1_ROI_definition_windows.csv'));

    save(fullfile(cfg.outFolder, 'E1_ROI_definition.mat'), ...
        'ROI', 'D', 'cfg', '-v7.3');

end

function EEG = reorder_eeg_to_template(EEG, templateLabels, templateTimes)

    labels = string({EEG.chanlocs.labels});
    idx = labels_to_indices(labels, templateLabels);

    if numel(idx) ~= numel(templateLabels)
        missing = templateLabels(~ismember(lower(templateLabels), lower(labels)));
        error('EEG is missing template channels: %s', strjoin(cellstr(missing), ', '));
    end

    if numel(EEG.times) ~= numel(templateTimes) || ...
            any(abs(double(EEG.times(:))' - double(templateTimes(:))') > 1e-6)
        error('EEG time vector does not match the E1 template.');
    end

    EEG.data = EEG.data(idx,:,:);
    EEG.chanlocs = EEG.chanlocs(idx);
    EEG.nbchan = numel(idx);

end

function idx = labels_to_indices(allLabels, requestedLabels)

    allLabels = string(allLabels(:));
    requestedLabels = string(requestedLabels(:));
    idx = zeros(0,1);

    for i = 1:numel(requestedLabels)
        hit = find(strcmpi(allLabels, requestedLabels(i)), 1, 'first');
        if ~isempty(hit)
            idx(end+1,1) = hit; %#ok<AGROW>
        end
    end

end

function [Lstd, Ldev, Gstd, Gdev] = compute_subject_lg_erps(EEG, targets)

    if ~(isfield(EEG, 'epoch') && isfield(EEG.epoch, 'eventlg_type'))
        error('EEG.epoch.eventlg_type is missing.');
    end

    M = false(4, EEG.trials);

    for ep = 1:EEG.trials
        labelsThisEpoch = to_string_list(EEG.epoch(ep).eventlg_type);
        for k = 1:4
            M(k,ep) = any(strcmpi(strtrim(labelsThisEpoch), targets(k)));
        end
    end

    idxFLS = find(M(1,:));
    idxRLS = find(M(2,:));
    idxFLD = find(M(3,:));
    idxRLD = find(M(4,:));

    if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
        error('At least one EEG local-global condition is empty.');
    end

    ERP_FLS = mean(double(EEG.data(:,:,idxFLS)), 3, 'omitnan');
    ERP_RLS = mean(double(EEG.data(:,:,idxRLS)), 3, 'omitnan');
    ERP_FLD = mean(double(EEG.data(:,:,idxFLD)), 3, 'omitnan');
    ERP_RLD = mean(double(EEG.data(:,:,idxRLD)), 3, 'omitnan');

    Lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3, 'omitnan');
    Ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3, 'omitnan');
    Gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3, 'omitnan');
    Gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3, 'omitnan');

end

function [pFWE, info] = fwe_maxT_paired_roi(devSub, stdSub, nPerm, seed, verbose, progressStep)
% Same paired sign-flip max-|t| FWE logic used in B7_EEG_stats.m.

    rng(seed, 'twister');

    D = double(devSub) - double(stdSub);
    nSub = size(D,1);

    tObs = paired_tmap_roi(D);
    absTobs = abs(tObs(:));

    maxT = zeros(nPerm,1);

    for ip = 1:nPerm

        sgn = ones(nSub,1);
        sgn(rand(nSub,1) > 0.5) = -1;

        Dp = D .* reshape(sgn, [nSub 1 1]);
        tPerm = paired_tmap_roi(Dp);
        maxT(ip) = max(abs(tPerm(:)));

        if verbose && mod(ip, progressStep) == 0
            fprintf('      FWE maxT permutation %d/%d\n', ip, nPerm);
        end

    end

    pVec = zeros(size(absTobs));
    for i = 1:numel(absTobs)
        pVec(i) = (1 + sum(maxT >= absTobs(i))) / (1 + nPerm);
    end

    pFWE = reshape(pVec, size(tObs));

    info = struct();
    info.N = nSub;
    info.df = nSub - 1;
    info.nPerm = nPerm;
    info.seed = seed;
    info.maxT = maxT;

end

function tMap = paired_tmap_roi(D)

    n = size(D,1);
    m = squeeze(mean(D,1));
    s = squeeze(std(D,0,1));
    tMap = m ./ (s ./ sqrt(n));
    tMap(~isfinite(tMap)) = 0;

end

function maskTime = triple_criterion_time_mask_roi(pMap, alpha, minElec, minConsec)

    sigElecCount = sum(pMap <= alpha, 1);
    initialMask = sigElecCount >= minElec;

    maskTime = false(size(initialMask));
    d = diff([false initialMask false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    for k = 1:numel(starts)
        if ends(k) - starts(k) + 1 >= minConsec
            maskTime(starts(k):ends(k)) = true;
        end
    end

end

function E = extract_erp_subject_features(eegFile, expKey, cfg, ROI) %#ok<INUSD>

    [filePath, fileName, fileExt] = fileparts(eegFile);

    EEG = pop_loadset('filename', [fileName fileExt], 'filepath', filePath);
    EEG = eeg_checkset(EEG);

    if ~(isfield(EEG, 'epoch') && isfield(EEG.epoch, 'eventlg_type'))
        error('EEG.epoch.eventlg_type is missing.');
    end

    labels = string({EEG.chanlocs.labels});

    localIdx = labels_to_indices(labels, ROI.Local.ChannelLabels);
    globalIdx = labels_to_indices(labels, ROI.Global.ChannelLabels);

    if isempty(localIdx)
        error('None of the selected Local ROI channels were found in %s.', eegFile);
    end
    if isempty(globalIdx)
        error('None of the selected Global ROI channels were found in %s.', eegFile);
    end

    M = false(4, EEG.trials);

    for ep = 1:EEG.trials

        labelsThisEpoch = to_string_list(EEG.epoch(ep).eventlg_type);

        for k = 1:4
            M(k,ep) = any(strcmpi(strtrim(labelsThisEpoch), cfg.targets(k)));
        end

    end

    idxFLS = find(M(1,:));
    idxRLS = find(M(2,:));
    idxFLD = find(M(3,:));
    idxRLD = find(M(4,:));

    if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
        error('At least one EEG local-global condition is empty.');
    end

    % Equal-weight elementary condition averages, matching B7_EEG_stats.m
    ERP_FLS = mean(double(EEG.data(:,:,idxFLS)), 3, 'omitnan');
    ERP_RLS = mean(double(EEG.data(:,:,idxRLS)), 3, 'omitnan');
    ERP_FLD = mean(double(EEG.data(:,:,idxFLD)), 3, 'omitnan');
    ERP_RLD = mean(double(EEG.data(:,:,idxRLD)), 3, 'omitnan');

    Lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3, 'omitnan');
    Ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3, 'omitnan');

    Gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3, 'omitnan');
    Gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3, 'omitnan');

    times = double(EEG.times(:))';

    % Average the condition waveforms spatially across the E1-derived ROI.
    roiLS = mean(Lstd(localIdx,:), 1, 'omitnan');
    roiLD = mean(Ldev(localIdx,:), 1, 'omitnan');

    roiGS = mean(Gstd(globalIdx,:), 1, 'omitnan');
    roiGD = mean(Gdev(globalIdx,:), 1, 'omitnan');

    E.Local = effect_features_from_curves( ...
        roiLS, roiLD, times, ROI.Local.WindowMS, cfg.peakDefinition);

    E.Global = effect_features_from_curves( ...
        roiGS, roiGD, times, ROI.Global.WindowMS, cfg.peakDefinition);

    E.Local.ChannelLabels  = ROI.Local.ChannelLabels;
    E.Global.ChannelLabels = ROI.Global.ChannelLabels;

end

function S = to_string_list(v)

    if isstring(v)
        S = v(:);
    elseif ischar(v)
        S = string(v);
    elseif iscell(v)
        S = strings(0,1);

        for i = 1:numel(v)
            if isstring(v{i}) || ischar(v{i})
                S(end+1,1) = string(v{i}); %#ok<AGROW>
            elseif isnumeric(v{i}) && isscalar(v{i})
                S(end+1,1) = string(v{i}); %#ok<AGROW>
            end
        end

    elseif isnumeric(v)
        S = string(v(:));
    else
        S = strings(0,1);
    end

    S = strtrim(regexprep(S, '\s+', ' '));

end

function F = effect_features_from_curves(stdCurve, devCurve, times, windowMS, peakDefinition)

    stdCurve = double(stdCurve(:))';
    devCurve = double(devCurve(:))';
    times    = double(times(:))';

    if numel(stdCurve) ~= numel(devCurve) || numel(stdCurve) ~= numel(times)
        error('Condition curves and time vector have incompatible lengths.');
    end

    idx = times >= windowMS(1) & times <= windowMS(2);

    if ~any(idx)
        error('No samples were found in window %.1f to %.1f ms.', ...
            windowMS(1), windowMS(2));
    end

    diffCurve = devCurve - stdCurve;

    F.Window = windowMS;
    F.StandardMean = mean(stdCurve(idx), 'omitnan');
    F.DeviantMean  = mean(devCurve(idx), 'omitnan');
    F.MeanEffect   = mean(diffCurve(idx), 'omitnan');

    wTimes = times(idx);
    wDiff  = diffCurve(idx);

    [F.PeakEffect, F.PeakLatency] = ...
        extract_peak(wDiff, wTimes, peakDefinition);

end

function [peakValue, peakLatency] = extract_peak(values, times, method)

    values = double(values(:));
    times  = double(times(:));

    valid = isfinite(values) & isfinite(times);
    values = values(valid);
    times  = times(valid);

    if isempty(values)
        peakValue = NaN;
        peakLatency = NaN;
        return;
    end

    switch lower(string(method))

        case "signed_abs_max"
            [~, idx] = max(abs(values));

        case "positive_max"
            [~, idx] = max(values);

        otherwise
            error('Unknown peak definition: %s', method);

    end

    peakValue = values(idx);
    peakLatency = times(idx);

end

function FEATURES = add_general_composites(FEATURES, mode)

    switch lower(string(mode))

        case "raw_equal_weight"

            FEATURES.PupilPeakGeneral = row_mean_two( ...
                FEATURES.PupilPeakLocal, FEATURES.PupilPeakGlobal);

            FEATURES.PupilMeanGeneral = row_mean_two( ...
                FEATURES.PupilMeanLocal, FEATURES.PupilMeanGlobal);

            FEATURES.ERPPeakGeneral = row_mean_two( ...
                FEATURES.ERPPeakLocal, FEATURES.ERPPeakGlobal);

            FEATURES.ERPMeanGeneral = row_mean_two( ...
                FEATURES.ERPMeanLocal, FEATURES.ERPMeanGlobal);

        case "zscore_equal_weight"

            FEATURES.PupilPeakGeneral = nan(height(FEATURES),1);
            FEATURES.PupilMeanGeneral = nan(height(FEATURES),1);
            FEATURES.ERPPeakGeneral   = nan(height(FEATURES),1);
            FEATURES.ERPMeanGeneral   = nan(height(FEATURES),1);

            exps = unique(FEATURES.Experiment);

            for e = 1:numel(exps)

                idx = FEATURES.Experiment == exps(e);

                FEATURES.PupilPeakGeneral(idx) = row_mean_two( ...
                    zscore_safe(FEATURES.PupilPeakLocal(idx)), ...
                    zscore_safe(FEATURES.PupilPeakGlobal(idx)));

                FEATURES.PupilMeanGeneral(idx) = row_mean_two( ...
                    zscore_safe(FEATURES.PupilMeanLocal(idx)), ...
                    zscore_safe(FEATURES.PupilMeanGlobal(idx)));

                FEATURES.ERPPeakGeneral(idx) = row_mean_two( ...
                    zscore_safe(FEATURES.ERPPeakLocal(idx)), ...
                    zscore_safe(FEATURES.ERPPeakGlobal(idx)));

                FEATURES.ERPMeanGeneral(idx) = row_mean_two( ...
                    zscore_safe(FEATURES.ERPMeanLocal(idx)), ...
                    zscore_safe(FEATURES.ERPMeanGlobal(idx)));

            end

        otherwise
            error('Unknown generalCompositeMode: %s', mode);

    end

end

function z = zscore_safe(x)

    x = double(x(:));
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');

    if ~isfinite(sd) || sd == 0
        z = nan(size(x));
    else
        z = (x - mu) ./ sd;
    end

end

function y = row_mean_two(a, b)

    a = double(a(:));
    b = double(b(:));

    y = nan(size(a));

    for i = 1:numel(a)
        vals = [a(i), b(i)];
        vals = vals(isfinite(vals));

        if ~isempty(vals)
            y(i) = mean(vals);
        end
    end

end

function ANALYSES = build_analysis_definitions()

    rows = struct([]);
    r = 0;

    expList = ["E1","E2"];
    contrastList = ["Local","Global","General"];

    % Paper-closest: pupil peak versus ERP mean
    for e = 1:numel(expList)
        for c = 1:numel(contrastList)

            r = r + 1;
            rows(r).AnalysisID = string(sprintf('%s_Alamia_%s', char(expList(e)), char(contrastList(c))));
            rows(r).AnalysisFamily = "Alamia";
            rows(r).Experiment = expList(e);
            rows(r).Measure = "PupilPeak_ERPMean";
            rows(r).Contrast = contrastList(c);
            rows(r).PupilVariable = "PupilPeak" + contrastList(c);
            rows(r).ERPVariable = "ERPMean" + contrastList(c);

        end
    end

    % Requested peak versus peak
    for e = 1:numel(expList)
        for c = 1:numel(contrastList)

            r = r + 1;
            rows(r).AnalysisID = string(sprintf('%s_PeakPeak_%s', char(expList(e)), char(contrastList(c))));
            rows(r).AnalysisFamily = "PeakPeak";
            rows(r).Experiment = expList(e);
            rows(r).Measure = "PeakPeak";
            rows(r).Contrast = contrastList(c);
            rows(r).PupilVariable = "PupilPeak" + contrastList(c);
            rows(r).ERPVariable = "ERPPeak" + contrastList(c);

        end
    end

    % Requested mean versus mean
    for e = 1:numel(expList)
        for c = 1:numel(contrastList)

            r = r + 1;
            rows(r).AnalysisID = string(sprintf('%s_MeanMean_%s', char(expList(e)), char(contrastList(c))));
            rows(r).AnalysisFamily = "MeanMean";
            rows(r).Experiment = expList(e);
            rows(r).Measure = "MeanMean";
            rows(r).Contrast = contrastList(c);
            rows(r).PupilVariable = "PupilMean" + contrastList(c);
            rows(r).ERPVariable = "ERPMean" + contrastList(c);

        end
    end

    ANALYSES = struct2table(rows);

end

function ST = run_correlation(x, y, subjects, cfg, analysisIndex)

    x = double(x(:));
    y = double(y(:));
    subjects = string(subjects(:));

    valid = isfinite(x) & isfinite(y);

    x = x(valid);
    y = y(valid);
    subjects = subjects(valid);

    ST = struct();
    ST.x = x;
    ST.y = y;
    ST.subjects = subjects;
    ST.N = numel(x);

    ST.rPearson = NaN;
    ST.pPearson = NaN;
    ST.pPermutation = NaN;
    ST.rSpearman = NaN;
    ST.pSpearman = NaN;
    ST.CI95_Pearson = [NaN NaN];

    if ST.N < 4 || numel(unique(x)) < 2 || numel(unique(y)) < 2
        return;
    end

    [ST.rPearson, ST.pPearson] = corr( ...
        x, y, 'Type', 'Pearson', 'Rows', 'complete');

    [ST.rSpearman, ST.pSpearman] = corr( ...
        x, y, 'Type', 'Spearman', 'Rows', 'complete');

    ST.CI95_Pearson = pearson_ci(ST.rPearson, ST.N, 0.05);

    rng(cfg.permutationSeed + analysisIndex, 'twister');

    rNull = nan(cfg.nPerm,1);

    for p = 1:cfg.nPerm
        yp = y(randperm(ST.N));
        rNull(p) = corr(x, yp, 'Type', 'Pearson', 'Rows', 'complete');
    end

    ST.pPermutation = ...
        (sum(abs(rNull) >= abs(ST.rPearson)) + 1) / (cfg.nPerm + 1);

end

function ci = pearson_ci(r, n, alpha)

    if n <= 3 || ~isfinite(r) || abs(r) >= 1
        ci = [NaN NaN];
        return;
    end

    z = atanh(r);
    se = 1 / sqrt(n - 3);
    zcrit = sqrt(2) * erfinv(1 - alpha);

    ci = tanh([z - zcrit*se, z + zcrit*se]);

end

function [BAYES_RESULTS, success, message] = ...
        run_bayesian_correlations_batch(BAYES_LONG, cfg)

    BAYES_RESULTS = table();
    success = false;
    message = "Bayesian analysis was not attempted.";

    if isempty(BAYES_LONG)
        message = "No complete correlation data were available.";
        return;
    end

    tempFolder = [tempname '_PupilERP_Bayes'];
    mkdir(tempFolder);

    inputCsv  = fullfile(tempFolder, 'bayes_input.csv');
    outputCsv = fullfile(tempFolder, 'bayes_output.csv');
    rFile     = fullfile(tempFolder, 'run_bayes_correlations.R');

    writetable(BAYES_LONG, inputCsv);

    fid = fopen(rFile, 'w');

    if fid < 0
        message = "Could not create the temporary R script.";
        return;
    end

    cleanupObj = onCleanup(@() safe_cleanup_temp(tempFolder)); %#ok<NASGU>

    fprintf(fid, '%s\n', ...
        'args <- commandArgs(trailingOnly=TRUE)', ...
        'input_file <- args[1]', ...
        'output_file <- args[2]', ...
        'iterations <- as.integer(args[3])', ...
        'if (!requireNamespace("BayesFactor", quietly=TRUE)) {', ...
        '  stop("The R package BayesFactor is not installed.")', ...
        '}', ...
        'set.seed(2026)', ...
        'd <- read.csv(input_file, stringsAsFactors=FALSE)', ...
        'ids <- unique(d$AnalysisID)', ...
        'out <- vector("list", length(ids))', ...
        'for (i in seq_along(ids)) {', ...
        '  id <- ids[i]', ...
        '  z <- d[d$AnalysisID == id, , drop=FALSE]', ...
        '  z <- z[is.finite(z$Pupil) & is.finite(z$ERP), , drop=FALSE]', ...
        '  bf_value <- NA_real_', ...
        '  ci_low <- NA_real_', ...
        '  ci_high <- NA_real_', ...
        '  err <- ""', ...
        '  if (nrow(z) >= 4 && length(unique(z$Pupil)) > 1 && length(unique(z$ERP)) > 1) {', ...
        '    tryCatch({', ...
        '      bf_obj <- BayesFactor::correlationBF(', ...
        '        y=z$ERP, x=z$Pupil, rscale="ultrawide")', ...
        '      bf_value <- as.numeric(BayesFactor::extractBF(', ...
        '        bf_obj, onlybf=TRUE))', ...
        '      post <- BayesFactor::posterior(', ...
        '        bf_obj, iterations=iterations, progress=FALSE)', ...
        '      rho <- as.numeric(post[, "rho"])', ...
        '      qs <- stats::quantile(rho, probs=c(0.025,0.975), na.rm=TRUE)', ...
        '      ci_low <- unname(qs[1])', ...
        '      ci_high <- unname(qs[2])', ...
        '    }, error=function(e) {', ...
        '      err <<- conditionMessage(e)', ...
        '    })', ...
        '  }', ...
        '  out[[i]] <- data.frame(', ...
        '    AnalysisID=id, BF10=bf_value,', ...
        '    CrI95_Low=ci_low, CrI95_High=ci_high,', ...
        '    Method="Jeffreys uniform prior via R BayesFactor",', ...
        '    Error=err, stringsAsFactors=FALSE)', ...
        '}', ...
        'out <- do.call(rbind, out)', ...
        'write.csv(out, output_file, row.names=FALSE)');

    fclose(fid);

    if isempty(cfg.rscriptPath)
        rscript = 'Rscript';
    else
        rscript = cfg.rscriptPath;
    end

    cmd = sprintf('"%s" "%s" "%s" "%s" %d', ...
        rscript, rFile, inputCsv, outputCsv, cfg.bayesPosteriorIterations);

    [status, cmdOut] = system(cmd);

    if status ~= 0
        message = "Exact R BayesFactor analysis failed: " + string(strtrim(cmdOut));
        return;
    end

    if exist(outputCsv, 'file') ~= 2
        message = "R finished without creating the Bayesian output CSV.";
        return;
    end

    try
        BAYES_RESULTS = readtable(outputCsv, 'TextType', 'string');
    catch
        message = "Could not read the Bayesian output.";
        return;
    end

    if ~all(ismember({'AnalysisID','BF10','CrI95_Low','CrI95_High','Method'}, ...
            BAYES_RESULTS.Properties.VariableNames))
        message = "The Bayesian output table did not contain the expected columns.";
        return;
    end

    if all(~isfinite(BAYES_RESULTS.BF10))
        message = "R ran, but no finite Bayes factors were returned.";
        return;
    end

    success = true;
    message = "Exact Jeffreys Bayesian correlation BF10 values were computed with R BayesFactor.";

end

function safe_cleanup_temp(folder)

    try
        if exist(folder, 'dir')
            rmdir(folder, 's');
        end
    catch
        % Do not interrupt the main analysis because temporary cleanup failed.
    end

end

function BF10 = bic_bf10_from_correlation(r, n)

    if ~isfinite(r) || ~isfinite(n) || n < 4
        BF10 = NaN;
        return;
    end

    r2 = min(max(r.^2, 0), 1 - eps);

    % BIC approximation for regression with one predictor:
    % BF10 approximately exp((BIC0 - BIC1)/2)
    logBF10 = 0.5 * (-n * log(1 - r2) - log(n));

    if logBF10 > log(realmax)
        BF10 = Inf;
    else
        BF10 = exp(logBF10);
    end

end

function pAdj = holm_adjust(p)

    p = double(p(:));
    pAdj = nan(size(p));

    valid = isfinite(p);
    pv = p(valid);

    if isempty(pv)
        return;
    end

    [pSorted, order] = sort(pv);
    m = numel(pSorted);

    adjSorted = (m - (1:m)' + 1) .* pSorted;
    adjSorted = cummax(adjSorted);
    adjSorted(adjSorted > 1) = 1;

    temp = nan(size(pv));
    temp(order) = adjSorted;

    pAdj(valid) = temp;

end

function fig = plot_local_global_figure(FEATURES, RESULTS, family, mainTitle, cfg)

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [80 80 1450 1000]);

    tl = tiledlayout(fig, 2, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    expList = ["E1","E2"];
    contrastList = ["Local","Global"];

    for e = 1:2
        for c = 1:2

            ax = nexttile(tl);

            analysisID = expList(e) + "_" + family + "_" + contrastList(c);

            plot_one_correlation_panel( ...
                ax, FEATURES, RESULTS, analysisID, cfg);

        end
    end

    sgtitle(tl, mainTitle, 'FontSize', 20, 'FontWeight', 'bold');

end

function fig = plot_general_figure(FEATURES, RESULTS, family, mainTitle, cfg)

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 120 1450 600]);

    tl = tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    expList = ["E1","E2"];

    for e = 1:2

        ax = nexttile(tl);

        analysisID = expList(e) + "_" + family + "_General";

        plot_one_correlation_panel( ...
            ax, FEATURES, RESULTS, analysisID, cfg);

    end

    sgtitle(tl, mainTitle, 'FontSize', 20, 'FontWeight', 'bold');

end

function plot_one_correlation_panel(ax, FEATURES, RESULTS, analysisID, cfg)

    ridx = find(RESULTS.AnalysisID == analysisID, 1, 'first');

    if isempty(ridx)
        title(ax, analysisID + " not found");
        axis(ax, 'off');
        return;
    end

    expKey = RESULTS.Experiment(ridx);
    xVar = RESULTS.PupilVariable(ridx);
    yVar = RESULTS.ERPVariable(ridx);
    contrast = RESULTS.Contrast(ridx);

    idx = FEATURES.Experiment == expKey;

    x = double(FEATURES.(char(xVar))(idx));
    y = double(FEATURES.(char(yVar))(idx));
    subj = FEATURES.Subject(idx);

    valid = isfinite(x) & isfinite(y);

    x = x(valid);
    y = y(valid);
    subj = subj(valid);

    hold(ax, 'on');

    scatter(ax, x, y, 72, ...
        'filled', ...
        'MarkerFaceColor', cfg.colors(1,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.82, ...
        'LineWidth', 0.75);

    if numel(x) >= 2 && numel(unique(x)) > 1

        b = polyfit(x, y, 1);
        xx = linspace(min(x), max(x), 200);
        yy = polyval(b, xx);

        plot(ax, xx, yy, ...
            'Color', cfg.colors(2,:), ...
            'LineWidth', 2.8);

    end

    if cfg.showSubjectLabels
        for i = 1:numel(x)
            text(ax, x(i), y(i), "  " + subj(i), ...
                'FontSize', 8, ...
                'Color', [0.25 0.25 0.25], ...
                'VerticalAlignment', 'middle');
        end
    end

    xlabel(ax, make_axis_label(xVar), ...
        'FontSize', 12, 'FontWeight', 'bold');

    ylabel(ax, make_axis_label(yVar), ...
        'FontSize', 12, 'FontWeight', 'bold');

    if expKey == "E1"
        expLabel = "Counting (E1)";
    else
        expLabel = "Visual interference (E2)";
    end

    title(ax, expLabel + " - " + contrast, ...
        'FontSize', 15, 'FontWeight', 'bold');

    statText = compose_stat_text(RESULTS(ridx,:));

    xl = padded_limits(x, 0.12);
    yl = padded_limits(y, 0.16);

    xlim(ax, xl);
    ylim(ax, yl);

    text(ax, ...
        xl(1) + 0.045 * range(xl), ...
        yl(2) - 0.055 * range(yl), ...
        statText, ...
        'FontSize', 10.5, ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.78 0.78 0.78], ...
        'Margin', 6);

    grid(ax, 'on');
    box(ax, 'off');

    set(ax, ...
        'FontSize', 11, ...
        'LineWidth', 1.2, ...
        'TickDir', 'out', ...
        'Layer', 'top');

    hold(ax, 'off');

end

function txt = compose_stat_text(resultRow)

    BF = resultRow.BF10;
    method = resultRow.BFmethod;

    if contains(lower(method), 'bic')
        bfLabel = "BF_{10} \approx " + format_bf(BF);
    else
        bfLabel = "BF_{10} = " + format_bf(BF);
    end

    txt = sprintf([ ...
        'N = %d\n' ...
        'r = %.2f\n' ...
        'p = %.3f\n' ...
        'perm p = %.3f\n' ...
        'p_{Holm} = %.3f\n' ...
        '%s'], ...
        resultRow.N, ...
        resultRow.rPearson, ...
        resultRow.pPearson, ...
        resultRow.pPermutation, ...
        resultRow.pPearsonHolm, ...
        char(bfLabel));

end

function lbl = make_axis_label(varName)

    v = string(varName);

    if startsWith(v, "PupilPeak")
        measure = "peak";
        unit = "z-score";
        contrast = extractAfter(v, "PupilPeak");
        lbl = sprintf('Pupil %s %s effect (%s)', char(lower(contrast)), measure, unit);

    elseif startsWith(v, "PupilMean")
        measure = "mean";
        unit = "z-score";
        contrast = extractAfter(v, "PupilMean");
        lbl = sprintf('Pupil %s %s effect (%s)', char(lower(contrast)), measure, unit);

    elseif startsWith(v, "ERPPeak")
        measure = "peak";
        unit = '\muV';
        contrast = extractAfter(v, "ERPPeak");
        lbl = sprintf('ROI ERP %s %s effect (%s)', char(lower(contrast)), measure, unit);

    elseif startsWith(v, "ERPMean")
        measure = "mean";
        unit = '\muV';
        contrast = extractAfter(v, "ERPMean");
        lbl = sprintf('ROI ERP %s %s effect (%s)', char(lower(contrast)), measure, unit);

    else
        lbl = char(replace(v, "_", " "));

    end

end

function lim = padded_limits(x, fraction)

    x = double(x(:));
    x = x(isfinite(x));

    if isempty(x)
        lim = [0 1];
        return;
    end

    xmin = min(x);
    xmax = max(x);

    if xmin == xmax
        pad = max(1, abs(xmin) * 0.2);
        lim = [xmin-pad, xmax+pad];
    else
        pad = fraction * (xmax - xmin);
        lim = [xmin-pad, xmax+pad];
    end

end

function save_figure_pair(fig, outFolder, baseName, cfg)

    if ~cfg.saveFigures
        return;
    end

    pngFile = fullfile(outFolder, [baseName '.png']);
    figFile = fullfile(outFolder, [baseName '.fig']);

    exportgraphics(fig, pngFile, 'Resolution', cfg.exportDPI);
    savefig(fig, figFile);

    fprintf('Saved figure:\n  %s\n', pngFile);

end

function s = format_bf(BF)

    if ~isfinite(BF)
        if isinf(BF)
            s = "Inf";
        else
            s = "NA";
        end
        return;
    end

    if BF >= 1000 || (BF > 0 && BF < 0.001)
        s = string(sprintf('%.2e', BF));
    elseif BF >= 100
        s = string(sprintf('%.1f', BF));
    elseif BF >= 10
        s = string(sprintf('%.2f', BF));
    else
        s = string(sprintf('%.3f', BF));
    end

end
