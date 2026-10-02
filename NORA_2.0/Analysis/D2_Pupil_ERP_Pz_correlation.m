clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% D2_Pupil_ERP_Pz_correlation.m
%
% PURPOSE
%   Participant-level correlations between pupil deviance effects and
%   Pz ERP deviance effects in the local-global paradigm.
%
%   The script loads the same subject-level preprocessed files used by:
%       C3_Pupil_stats.m
%       B7_EEG_stats.m
%
%   It computes, for each participant and experiment:
%       LOCAL pupil effect  = Local Deviant  - Local Standard
%       GLOBAL pupil effect = Global Deviant - Global Standard
%       LOCAL Pz ERP effect = Local Deviant  - Local Standard
%       GLOBAL Pz ERP effect= Global Deviant - Global Standard
%
% PRIMARY / PAPER-CLOSEST ANALYSIS
%   Alamia et al. Figure 4D did NOT correlate an ERP peak with a pupil peak.
%   They correlated:
%       - the signed absolute pupil maximum representing the rare-minus-
%         frequent pupil response
%       - the MEAN rare-minus-frequent ERP difference, averaged over the
%         second significant EEG cluster and its time interval.
%
%   Therefore this script creates:
%
%   A) TWO ALAMIA-LIKE FIGURES
%      1. Pupil PEAK effect vs Pz ERP MEAN effect:
%         local/global panels, separately for E1 and E2
%      2. The same for the general local-global composite
%
%   B) FOUR REQUESTED EXPLORATORY FIGURES
%      3. Pupil PEAK vs Pz ERP PEAK: local/global
%      4. Pupil PEAK vs Pz ERP PEAK: general composite
%      5. Pupil MEAN vs Pz ERP MEAN: local/global
%      6. Pupil MEAN vs Pz ERP MEAN: general composite
%
%   Set cfg.makeAlamiaLikeFigures = false if you want only the four
%   requested matched-measure figures.
%
% TIME WINDOWS
%   Pupil windows are the significant cluster windows from the manuscript:
%       E1 Local  = 44-1232 ms
%       E1 Global = 288-2000 ms
%       E2 Local  = 0-690 ms
%       E2 Global = 446-1818 ms
%
%   Pz ERP windows are the significant E1 counting-derived windows:
%       Local  = 236-308 ms
%       Global = 208-532 ms
%
%   As requested, the SAME E1-derived ERP windows are applied to E2,
%   because E2 had no corrected group-level ERP effect.
%
% BAYESIAN CORRELATION
%   The primary coefficient is Pearson's r, matching Alamia et al.
%
%   The script first tries to obtain an exact Jeffreys Bayesian correlation
%   BF10 using the R package BayesFactor with an ultrawide rscale of 1.
%   This corresponds to a uniform prior on rho from -1 to 1, close to the
%   default JASP analysis used in Alamia et al.
%
%   If Rscript or the R package BayesFactor is unavailable, the script
%   automatically uses a clearly labelled BIC approximation to BF10.
%
% OUTPUT
%   - Subject-level extracted feature table (.csv)
%   - Full correlation-results table (.csv)
%   - MATLAB results file (.mat)
%   - Six figures by default (.png and .fig)
%
% IMPORTANT INTERPRETATION
%   The "General" effect is an exploratory equal-weight composite:
%       General effect = 0.5 * Local effect + 0.5 * Global effect
%
%   It is NOT a separate factorial contrast in the local-global design.
%   The local and global correlations should remain the main analyses.
%
% Miguel Contreras-Altamirano, 2026

%% ===================== USER SETTINGS =====================

cfg = struct();

% Main data paths
cfg.eeglabPath = 'L:\Downloads\eeglab2023.0\';
cfg.outroot    = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';
cfg.outFolder  = fullfile(cfg.outroot, 'GRAND_AVG', 'Pupil_ERP_correlation');

if ~exist(cfg.outFolder, 'dir')
    mkdir(cfg.outFolder);
end

% Data definitions
cfg.pupilField = 'PUPIL_Z';
cfg.erpChannel = 'Pz';

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

% E1 counting-derived Pz ERP windows.
% These windows are intentionally reused for E2.
cfg.erpWindows.Local  = [236 308];
cfg.erpWindows.Global = [208 532];

% Peak definition for difference waves:
%   'signed_abs_max' = sample with the largest absolute magnitude, retaining
%                      its sign. This is closest to Alamia's pupil method.
%   'positive_max'   = largest positive sample in the window.
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

% Number of posterior samples for the Bayesian 95% credible interval
cfg.bayesPosteriorIterations = 30000;

% If true, stop rather than using the BIC approximation when exact R-based
% BayesFactor analysis is unavailable.
cfg.requireExactBayesFactor = false;

% Figures
cfg.makeAlamiaLikeFigures = true;
cfg.saveFigures           = true;
cfg.showSubjectLabels     = false;
cfg.exportDPI             = 300;

% Requested colors
cfg.colors = [ ...
    0.00 0.45 0.85;   % participant dots: blue
    0.90 0.45 0.10];  % regression line: orange

% General composite
% 'raw_equal_weight': mean of raw Local and Global effect features.
% 'zscore_equal_weight': mean of within-experiment standardized Local and
%                        Global effects. The raw version is closer to
%                        "all standards versus all deviants".
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

%% ===================== EXTRACT SUBJECT FEATURES =====================

FEATURES = table();
expList = ["E1","E2"];

for e = 1:numel(expList)

    expKey = expList(e);

    fprintf('\n============================================================\n');
    fprintf('EXTRACTING PUPIL AND Pz ERP FEATURES: %s\n', expKey);
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
            E = extract_erp_subject_features(eegFile, expKey, cfg);
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

featureCsv = fullfile(cfg.outFolder, 'Pupil_ERP_subject_features.csv');
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

resultsCsv = fullfile(cfg.outFolder, 'Pupil_ERP_correlation_results.csv');
writetable(RESULTS, resultsCsv);

resultsMat = fullfile(cfg.outFolder, 'Pupil_ERP_correlation_results.mat');
save(resultsMat, 'FEATURES', 'RESULTS', 'ANALYSES', 'cfg', '-v7.3');

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
        'Pupil peak effect vs Pz ERP mean effect', ...
        cfg);

    save_figure_pair(fig1, cfg.outFolder, ...
        'Fig_PupilPeak_vs_ERPMean_LOCAL_GLOBAL_AlamiaLike', cfg);

    % 2. Paper-closest: pupil peak versus ERP mean, general composite
    fig2 = plot_general_figure( ...
        FEATURES, RESULTS, ...
        "Alamia", ...
        'Pupil peak effect vs Pz ERP mean effect: general composite', ...
        cfg);

    save_figure_pair(fig2, cfg.outFolder, ...
        'Fig_PupilPeak_vs_ERPMean_GENERAL_AlamiaLike', cfg);

end

% 3. Requested peak-peak local/global
fig3 = plot_local_global_figure( ...
    FEATURES, RESULTS, ...
    "PeakPeak", ...
    'Pupil peak effect vs Pz ERP peak effect', ...
    cfg);

save_figure_pair(fig3, cfg.outFolder, ...
    'Fig_PupilPeak_vs_ERPPeak_LOCAL_GLOBAL', cfg);

% 4. Requested peak-peak general
fig4 = plot_general_figure( ...
    FEATURES, RESULTS, ...
    "PeakPeak", ...
    'Pupil peak effect vs Pz ERP peak effect: general composite', ...
    cfg);

save_figure_pair(fig4, cfg.outFolder, ...
    'Fig_PupilPeak_vs_ERPPeak_GENERAL', cfg);

% 5. Requested mean-mean local/global
fig5 = plot_local_global_figure( ...
    FEATURES, RESULTS, ...
    "MeanMean", ...
    'Pupil mean effect vs Pz ERP mean effect', ...
    cfg);

save_figure_pair(fig5, cfg.outFolder, ...
    'Fig_PupilMean_vs_ERPMean_LOCAL_GLOBAL', cfg);

% 6. Requested mean-mean general
fig6 = plot_general_figure( ...
    FEATURES, RESULTS, ...
    "MeanMean", ...
    'Pupil mean effect vs Pz ERP mean effect: general composite', ...
    cfg);

save_figure_pair(fig6, cfg.outFolder, ...
    'Fig_PupilMean_vs_ERPMean_GENERAL', cfg);

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

function E = extract_erp_subject_features(eegFile, expKey, cfg) %#ok<INUSD>

    [filePath, fileName, fileExt] = fileparts(eegFile);

    EEG = pop_loadset('filename', [fileName fileExt], 'filepath', filePath);
    EEG = eeg_checkset(EEG);

    labels = string({EEG.chanlocs.labels});
    chIdx = find(strcmpi(labels, cfg.erpChannel), 1, 'first');

    if isempty(chIdx)
        error('ERP channel %s was not found.', cfg.erpChannel);
    end

    if ~(isfield(EEG, 'epoch') && isfield(EEG.epoch, 'eventlg_type'))
        error('EEG.epoch.eventlg_type is missing.');
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

    pzLS = squeeze(Lstd(chIdx,:))';
    pzLD = squeeze(Ldev(chIdx,:))';
    pzGS = squeeze(Gstd(chIdx,:))';
    pzGD = squeeze(Gdev(chIdx,:))';

    E.Local = effect_features_from_curves( ...
        pzLS, pzLD, times, cfg.erpWindows.Local, cfg.peakDefinition);

    E.Global = effect_features_from_curves( ...
        pzGS, pzGD, times, cfg.erpWindows.Global, cfg.peakDefinition);

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
        lbl = sprintf('Pz ERP %s %s effect (%s)', char(lower(contrast)), measure, unit);

    elseif startsWith(v, "ERPMean")
        measure = "mean";
        unit = '\muV';
        contrast = extractAfter(v, "ERPMean");
        lbl = sprintf('Pz ERP %s %s effect (%s)', char(lower(contrast)), measure, unit);

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
