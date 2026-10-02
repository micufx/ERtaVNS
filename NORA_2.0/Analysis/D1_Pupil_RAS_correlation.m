clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% C4_Pupil_RAS_GLOBAL_correlation.m

% PURPOSE
%   Participant-level correlation between GLOBAL pupil effect and GLOBAL
%   awareness score for BOTH experiments.
%
% MAIN ANALYSIS
%   E1: pupil global effect vs RAS_E1_Global
%   E2: pupil global effect vs RAS_E2_Global
%
% DEFINITIONS
%   Global pupil effect = Global Deviant - Global Standard
%                       = Rare/Frequent deviance effect
%
%   RAS_Global = Q3 + Q4
%   This is preferred over total RAS because Q1-Q2 mainly reflect more
%   general/local awareness and often show ceiling effects.
%
% INPUT
%   Pupil_RAS.mat
%   Expected to contain S.Combined with subject-level pupil condition values
%   and RAS variables.
%
% OUTPUT
%   - CORR_GlobalPupil_vs_GlobalRAS_Spearman.csv
%   - CORR_GlobalPupil_vs_GlobalRAS_Spearman.mat
%   - CORR_diagnostic_effects_and_RAS.csv
%   - Fig3E_style_E1_global_effect_vs_RAS_E1_Global.png/.fig
%   - Fig3E_style_E2_global_effect_vs_RAS_E2_Global.png/.fig
%   - Fig3E_style_BOTH_global_effect_vs_global_RAS.png/.fig
%
% IMPORTANT
%   The script FIRST tries to recompute the global pupil effect from the
%   condition means:
%       Global Deviant - Global Standard
%
%   If it cannot find the condition columns, it falls back to an existing
%   E1_global_effect / E2_global_effect column, but it prints a warning.
%
%   Check the command window after running to confirm which columns were
%   used for Global Standard and Global Deviant.

%% ===================== USER SETTINGS =====================

basepath  = 'L:\Downloads\scripts_Vagus\matlab\testing\out\GRAND_AVG\tables\Pupil_RAS';
matFile   = fullfile(basepath, 'Pupil_RAS.mat');
outFolder = fullfile(basepath, 'GRAND_AVG');

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% Correlation type:
% Spearman is recommended because RAS_Global is ordinal: 0, 1, or 2.
corrType = 'Spearman';      % 'Spearman' or 'Pearson'

% Permutation p-value for the correlation
doPermutationP = true;
nPerm = 10000;
rng(1, 'twister');

% If true, the script stops if it cannot recompute GD-GS from condition means.
% If false, it falls back to existing E1_global_effect / E2_global_effect.
requireConditionMeans = false;

% Figures
saveFigures = true;
showSubjectLabels = true;   % useful for checking outliers; set false for cleaner final figure
jitterY_forPlotOnly = true; % only visual jitter because RAS scores are discrete
jitterAmount = 0.045;

% Optional secondary check:
% This tests global pupil effect vs total RAS, but it is NOT the preferred
% main analysis because Total RAS mixes local/general and global items.
runTotalRASCheck = true;

%% ===================== LOAD DATA =====================

if ~exist(matFile, 'file')
    error('MAT file not found:\n%s', matFile);
end

S = load(matFile);
D = load_combined_as_table(S);
D = ensure_subject_column(D);

fprintf('\nLoaded table with %d rows and %d variables.\n', height(D), width(D));

if ismember('Has_RAS', D.Properties.VariableNames)
    fprintf('Subjects with matching RAS: %d / %d\n', sum(D.Has_RAS == 1), height(D));
end

%% ===================== ADD / STANDARDIZE RAS SCORES =====================

D = add_ras_scores(D);

% Main global awareness scores
D.RAS_E1_Global = D.Q3_E1 + D.Q4_E1;
D.RAS_E2_Global = D.Q3_E2 + D.Q4_E2;

% Total RAS scores computed directly from Q1-Q4
D.RAS_E1_Total_fromQ = D.Q1_E1 + D.Q2_E1 + D.Q3_E1 + D.Q4_E1;
D.RAS_E2_Total_fromQ = D.Q1_E2 + D.Q2_E2 + D.Q3_E2 + D.Q4_E2;

%% ===================== COMPUTE GLOBAL PUPIL EFFECTS =====================

% This recomputes:
%   E1_global_effect = E1_GlobalDeviant - E1_GlobalStandard
%   E2_global_effect = E2_GlobalDeviant - E2_GlobalStandard
[D, effectInfo] = add_global_pupil_effect(D, 'E1', requireConditionMeans);
[D, effectInfo] = add_global_pupil_effect(D, 'E2', requireConditionMeans, effectInfo);

fprintf('\nGlobal effect definition used in this script:\n');
fprintf('  E1_global_effect = %s - %s\n', effectInfo.E1_GD, effectInfo.E1_GS);
fprintf('  E2_global_effect = %s - %s\n', effectInfo.E2_GD, effectInfo.E2_GS);
fprintf('  Positive values mean larger pupil response for Global Deviants than Global Standards.\n\n');

%% ===================== SAVE DIAGNOSTIC TABLE =====================

diagVars = { ...
    'Subject', ...
    'Q1_E1','Q2_E1','Q3_E1','Q4_E1','RAS_E1_Global','RAS_E1_Total_fromQ', ...
    'Q1_E2','Q2_E2','Q3_E2','Q4_E2','RAS_E2_Global','RAS_E2_Total_fromQ', ...
    'E1_global_standard_mean','E1_global_deviant_mean','E1_global_effect', ...
    'E2_global_standard_mean','E2_global_deviant_mean','E2_global_effect'};

diagVars = diagVars(ismember(diagVars, D.Properties.VariableNames));
Tdiag = D(:, diagVars);

outDiag = fullfile(outFolder, 'CORR_diagnostic_effects_and_RAS.csv');
writetable(Tdiag, outDiag);
fprintf('Saved diagnostic table:\n  %s\n', outDiag);

%% ===================== MAIN ANALYSES: GLOBAL PUPIL EFFECT VS GLOBAL RAS =====================

Analysis = struct( ...
    'experiment', {'E1', 'E2'}, ...
    'xVar',       {'E1_global_effect', 'E2_global_effect'}, ...
    'yVar',       {'RAS_E1_Global', 'RAS_E2_Global'}, ...
    'titleStr',   {'Counting (E1)', 'Visual interference (E2)'});

R = table();

for ia = 1:numel(Analysis)

    x = D.(Analysis(ia).xVar);
    y = D.(Analysis(ia).yVar);
    subj = string(D.Subject);

    ST = run_corr_with_permutation(x, y, subj, corrType, doPermutationP, nPerm);

    R.Experiment(ia,1) = string(Analysis(ia).experiment);
    R.Xvar(ia,1) = string(Analysis(ia).xVar);
    R.Yvar(ia,1) = string(Analysis(ia).yVar);
    R.N(ia,1) = ST.N;
    R.rho_or_r(ia,1) = ST.r;
    R.p(ia,1) = ST.p;
    R.pPerm(ia,1) = ST.pPerm;
    R.CI95_low(ia,1) = ST.CI95(1);
    R.CI95_high(ia,1) = ST.CI95(2);
    R.N_unique_Y(ia,1) = numel(unique(ST.y));

    fprintf('\n%s: %s vs %s\n', Analysis(ia).experiment, Analysis(ia).xVar, Analysis(ia).yVar);
    fprintf('  N = %d\n', ST.N);
    fprintf('  %s = %.3f, p = %.4f, permutation p = %.4f\n', corrType, ST.r, ST.p, ST.pPerm);

end

% Holm correction across the two planned correlations: E1 and E2
R.pHolm = holm_adjust(R.p);
R.pPermHolm = holm_adjust(R.pPerm);

outCsv = fullfile(outFolder, sprintf('CORR_GlobalPupil_vs_GlobalRAS_%s.csv', corrType));
writetable(R, outCsv);
fprintf('\nSaved main correlation results:\n  %s\n', outCsv);

outMat = fullfile(outFolder, sprintf('CORR_GlobalPupil_vs_GlobalRAS_%s.mat', corrType));
save(outMat, 'R', 'D', 'Analysis', 'effectInfo');
fprintf('Saved main correlation MAT:\n  %s\n', outMat);

%% ===================== MAIN FIGURES: ONE FIGURE PER EXPERIMENT =====================

for ia = 1:numel(Analysis)

    x = D.(Analysis(ia).xVar);
    y = D.(Analysis(ia).yVar);
    subj = string(D.Subject);

    ST = run_corr_with_permutation(x, y, subj, corrType, doPermutationP, nPerm);
    ST.pHolm = R.pHolm(ia);
    ST.pPermHolm = R.pPermHolm(ia);

    fig = figure('Color','w','Position',[200 200 640 560]);
    ax = axes(fig);

    plot_corr_axis(ax, ST.x, ST.y, ST.subjects, ST, ...
        sprintf('Pupil global effect vs global awareness - %s', Analysis(ia).titleStr), ...
        make_pretty_label(Analysis(ia).xVar), ...
        make_pretty_label(Analysis(ia).yVar), ...
        corrType, showSubjectLabels, jitterY_forPlotOnly, jitterAmount);

    if saveFigures
        outPng = fullfile(outFolder, sprintf('Fig3E_style_%s_global_effect_vs_%s.png', ...
            Analysis(ia).experiment, Analysis(ia).yVar));
        outFig = fullfile(outFolder, sprintf('Fig3E_style_%s_global_effect_vs_%s.fig', ...
            Analysis(ia).experiment, Analysis(ia).yVar));

        exportgraphics(fig, outPng, 'Resolution', 300);
        savefig(fig, outFig);

        fprintf('Saved figure:\n  %s\n', outPng);
    end

end

%% ===================== COMBINED 2-PANEL FIGURE =====================

fig = figure('Color','w','Position',[100 100 1200 540]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for ia = 1:numel(Analysis)

    ax = nexttile;

    x = D.(Analysis(ia).xVar);
    y = D.(Analysis(ia).yVar);
    subj = string(D.Subject);

    ST = run_corr_with_permutation(x, y, subj, corrType, doPermutationP, nPerm);
    ST.pHolm = R.pHolm(ia);
    ST.pPermHolm = R.pPermHolm(ia);

    plot_corr_axis(ax, ST.x, ST.y, ST.subjects, ST, ...
        Analysis(ia).titleStr, ...
        make_pretty_label(Analysis(ia).xVar), ...
        make_pretty_label(Analysis(ia).yVar), ...
        corrType, showSubjectLabels, jitterY_forPlotOnly, jitterAmount);

end

sgtitle('Global pupil effect vs global awareness', 'FontSize', 17, 'FontWeight','bold');

if saveFigures
    outPng = fullfile(outFolder, 'Fig3E_style_BOTH_global_effect_vs_global_RAS.png');
    outFig = fullfile(outFolder, 'Fig3E_style_BOTH_global_effect_vs_global_RAS.fig');

    exportgraphics(fig, outPng, 'Resolution', 300);
    savefig(fig, outFig);

    fprintf('Saved combined figure:\n  %s\n', outPng);
end

%% ===================== OPTIONAL: TOTAL RAS CHECK =====================

if runTotalRASCheck

    fprintf('\nRunning optional secondary check: global pupil effect vs TOTAL RAS.\n');
    fprintf('Interpret this only as exploratory/descriptive because total RAS mixes local/general and global awareness.\n');

    TotalAnalysis = struct( ...
        'experiment', {'E1', 'E2'}, ...
        'xVar',       {'E1_global_effect', 'E2_global_effect'}, ...
        'yVar',       {'RAS_E1_Total_fromQ', 'RAS_E2_Total_fromQ'}, ...
        'titleStr',   {'Counting (E1)', 'Visual interference (E2)'});

    Rtotal = table();

    for ia = 1:numel(TotalAnalysis)

        x = D.(TotalAnalysis(ia).xVar);
        y = D.(TotalAnalysis(ia).yVar);
        subj = string(D.Subject);

        ST = run_corr_with_permutation(x, y, subj, corrType, doPermutationP, nPerm);

        Rtotal.Experiment(ia,1) = string(TotalAnalysis(ia).experiment);
        Rtotal.Xvar(ia,1) = string(TotalAnalysis(ia).xVar);
        Rtotal.Yvar(ia,1) = string(TotalAnalysis(ia).yVar);
        Rtotal.N(ia,1) = ST.N;
        Rtotal.rho_or_r(ia,1) = ST.r;
        Rtotal.p(ia,1) = ST.p;
        Rtotal.pPerm(ia,1) = ST.pPerm;
        Rtotal.CI95_low(ia,1) = ST.CI95(1);
        Rtotal.CI95_high(ia,1) = ST.CI95(2);
        Rtotal.N_unique_Y(ia,1) = numel(unique(ST.y));

    end

    Rtotal.pHolm = holm_adjust(Rtotal.p);
    Rtotal.pPermHolm = holm_adjust(Rtotal.pPerm);

    outTotalCsv = fullfile(outFolder, sprintf('CORR_GlobalPupil_vs_TotalRAS_%s_EXPLORATORY.csv', corrType));
    writetable(Rtotal, outTotalCsv);
    fprintf('Saved optional total-RAS correlation results:\n  %s\n', outTotalCsv);

end

fprintf('\nDone.\n');

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function D = load_combined_as_table(S)

    if isfield(S, 'Combined')
        A = S.Combined;
    else
        f = fieldnames(S);
        A = [];
        for i = 1:numel(f)
            if istable(S.(f{i})) || isstruct(S.(f{i}))
                A = S.(f{i});
                fprintf('No variable named Combined found. Using S.%s instead.\n', f{i});
                break;
            end
        end
        if isempty(A)
            error('No table or struct found inside the MAT file.');
        end
    end

    if istable(A)
        D = A;
    elseif isstruct(A)
        D = matstruct_to_table(A);
    else
        error('The loaded object is neither a table nor a struct.');
    end

end

function T = matstruct_to_table(A)

    f = fieldnames(A);
    T = table();

    for i = 1:numel(f)
        v = A.(f{i});

        if iscell(v)
            try
                T.(f{i}) = string(v(:));
            catch
                T.(f{i}) = v(:);
            end
        elseif ischar(v)
            T.(f{i}) = string(cellstr(v));
        elseif isstring(v)
            T.(f{i}) = v(:);
        elseif iscategorical(v)
            T.(f{i}) = string(v(:));
        elseif islogical(v)
            T.(f{i}) = double(v(:));
        else
            T.(f{i}) = double(v(:));
        end
    end

end

function D = ensure_subject_column(D)

    vars = D.Properties.VariableNames;

    if ismember('Subject', vars)
        subj = D.Subject;
    else
        idx = find(contains(lower(vars), 'subject') | contains(lower(vars), 'participant'), 1, 'first');
        if isempty(idx)
            warning('No Subject column found. Creating Subject labels by row order: S01-Sxx.');
            subj = (1:height(D))';
        else
            subj = D.(vars{idx});
            fprintf('Using %s as Subject column.\n', vars{idx});
        end
    end

    D.Subject = standardize_subject_labels(subj);

end

function subj = standardize_subject_labels(subjIn)

    if isnumeric(subjIn)
        subj = strings(numel(subjIn),1);
        for i = 1:numel(subjIn)
            subj(i) = sprintf('S%02d', subjIn(i));
        end
        return;
    end

    subjIn = string(subjIn);
    subj = strings(numel(subjIn),1);

    for i = 1:numel(subjIn)
        s = strtrim(subjIn(i));
        tok = regexp(char(s), '\d+', 'match');

        if ~isempty(tok)
            n = str2double(tok{end});
            if isfinite(n)
                subj(i) = sprintf('S%02d', n);
            else
                subj(i) = s;
            end
        else
            subj(i) = s;
        end
    end

end

function D = add_ras_scores(D)

    % Find and standardize Q1-Q4 variables for E1 and E2
    qNames = {'Q1','Q2','Q3','Q4'};
    expNames = {'E1','E2'};

    for ie = 1:numel(expNames)
        for iq = 1:numel(qNames)

            expTag = expNames{ie};
            qTag = qNames{iq};

            canonicalName = sprintf('%s_%s', qTag, expTag);

            v = find_ras_question_var(D, expTag, qTag);
            if isempty(v)
                error('Could not find RAS variable for %s %s.', expTag, qTag);
            end

            D.(canonicalName) = as_numeric(D.(v));

        end
    end

end

function v = find_ras_question_var(D, expTag, qTag)

    vars = D.Properties.VariableNames;

    aliases = { ...
        sprintf('%s_%s', qTag, expTag), ...
        sprintf('%s_%s', expTag, qTag), ...
        sprintf('%s%s', qTag, expTag), ...
        sprintf('%s%s', expTag, qTag)};

    v = find_var_exact_or_normalized(vars, aliases);
    if ~isempty(v)
        return;
    end

    expNorm = lower(expTag);
    qNorm = lower(qTag);

    bestScore = -Inf;
    bestVar = '';

    for i = 1:numel(vars)
        n = norm_name(vars{i});

        if contains(n, expNorm) && contains(n, qNorm)
            score = 1;
            if startsWith(n, qNorm), score = score + 2; end
            if startsWith(n, expNorm), score = score + 1; end

            if score > bestScore
                bestScore = score;
                bestVar = vars{i};
            end
        end
    end

    v = bestVar;

end

function [D, effectInfo] = add_global_pupil_effect(D, expTag, requireConditionMeans, effectInfo)

    if nargin < 4
        effectInfo = struct();
    end

    gsVar = find_best_global_condition_var(D, expTag, 'GS');
    gdVar = find_best_global_condition_var(D, expTag, 'GD');

    standardOut = sprintf('%s_global_standard_mean', expTag);
    deviantOut  = sprintf('%s_global_deviant_mean', expTag);
    effectOut   = sprintf('%s_global_effect', expTag);

    if ~isempty(gsVar) && ~isempty(gdVar)

        D.(standardOut) = as_numeric(D.(gsVar));
        D.(deviantOut)  = as_numeric(D.(gdVar));

        % Correct global effect:
        % Global Deviant - Global Standard
        D.(effectOut) = D.(deviantOut) - D.(standardOut);

        fprintf('\n%s global effect recomputed from condition means:\n', expTag);
        fprintf('  Global Standard = %s\n', gsVar);
        fprintf('  Global Deviant  = %s\n', gdVar);
        fprintf('  Effect          = %s - %s\n', gdVar, gsVar);

        effectInfo.(sprintf('%s_GS', expTag)) = gsVar;
        effectInfo.(sprintf('%s_GD', expTag)) = gdVar;
        effectInfo.(sprintf('%s_effect_source', expTag)) = 'recomputed_GD_minus_GS';

    else

        if requireConditionMeans
            fprintf('\nAvailable variable names:\n');
            disp(D.Properties.VariableNames');
            error(['Could not find Global Standard and Global Deviant condition columns for %s.\n' ...
                   'Set requireConditionMeans = false to fall back to existing %s_global_effect, or rename your columns.'], ...
                   expTag, expTag);
        end

        existingEffect = find_existing_global_effect_var(D, expTag);

        if isempty(existingEffect)
            fprintf('\nAvailable variable names:\n');
            disp(D.Properties.VariableNames');
            error('Could not recompute or find existing %s_global_effect.', expTag);
        end

        D.(effectOut) = as_numeric(D.(existingEffect));

        warning(['Could not find both condition means for %s. Falling back to existing column: %s.\n' ...
                 'Please confirm this existing column was computed as Global Deviant - Global Standard.'], ...
                 expTag, existingEffect);

        D.(standardOut) = NaN(height(D),1);
        D.(deviantOut) = NaN(height(D),1);

        effectInfo.(sprintf('%s_GS', expTag)) = 'not found';
        effectInfo.(sprintf('%s_GD', expTag)) = 'not found';
        effectInfo.(sprintf('%s_effect_source', expTag)) = existingEffect;

    end

end

function v = find_best_global_condition_var(D, expTag, condTag)

    vars = D.Properties.VariableNames;

    if strcmpi(condTag, 'GS')
        aliases = { ...
            sprintf('%s_GS_mean', expTag), ...
            sprintf('%s_GlobalStandard_mean', expTag), ...
            sprintf('%s_Global_Standard_mean', expTag), ...
            sprintf('%s_global_standard_mean', expTag), ...
            sprintf('%s_mean_GS', expTag), ...
            sprintf('%s_pupil_GS_mean', expTag)};
    else
        aliases = { ...
            sprintf('%s_GD_mean', expTag), ...
            sprintf('%s_GlobalDeviant_mean', expTag), ...
            sprintf('%s_Global_Deviant_mean', expTag), ...
            sprintf('%s_global_deviant_mean', expTag), ...
            sprintf('%s_mean_GD', expTag), ...
            sprintf('%s_pupil_GD_mean', expTag)};
    end

    v = find_var_exact_or_normalized(vars, aliases);
    if ~isempty(v)
        return;
    end

    expNorm = lower(expTag);
    bestScore = -Inf;
    bestVar = '';

    excludeTokens = {'ras','q1','q2','q3','q4','total','perceived','block', ...
                     'subject','has','original','effect','delta','latency', ...
                     'lat','peak','erp'};

    for i = 1:numel(vars)

        n = norm_name(vars{i});

        if ~contains(n, expNorm)
            continue;
        end

        if any_contains(n, excludeTokens)
            continue;
        end

        if strcmpi(condTag, 'GS')
            condMatch = contains(n, 'globalstandard') || ...
                        contains(n, [expNorm 'gs']) || ...
                        contains(n, 'gsmean') || ...
                        endsWith(n, 'gs');
        else
            condMatch = contains(n, 'globaldeviant') || ...
                        contains(n, [expNorm 'gd']) || ...
                        contains(n, 'gdmean') || ...
                        endsWith(n, 'gd');
        end

        if ~condMatch
            continue;
        end

        score = 0;

        if strcmpi(condTag, 'GS') && contains(n, 'globalstandard')
            score = score + 100;
        elseif strcmpi(condTag, 'GD') && contains(n, 'globaldeviant')
            score = score + 100;
        end

        if contains(n, 'mean'),  score = score + 20; end
        if contains(n, 'pupil'), score = score + 10; end
        if contains(n, 'z'),     score = score + 2;  end
        if startsWith(n, expNorm), score = score + 5; end

        if score > bestScore
            bestScore = score;
            bestVar = vars{i};
        end

    end

    v = bestVar;

end

function v = find_existing_global_effect_var(D, expTag)

    vars = D.Properties.VariableNames;

    aliases = { ...
        sprintf('%s_global_effect', expTag), ...
        sprintf('%s_GlobalEffect', expTag), ...
        sprintf('%s_global_pupil_effect', expTag), ...
        sprintf('%s_pupil_global_effect', expTag)};

    v = find_var_exact_or_normalized(vars, aliases);
    if ~isempty(v)
        return;
    end

    expNorm = lower(expTag);
    bestScore = -Inf;
    bestVar = '';

    for i = 1:numel(vars)
        n = norm_name(vars{i});

        if contains(n, expNorm) && contains(n, 'global') && contains(n, 'effect') && ...
                ~contains(n, 'delta') && ~contains(n, 'ras') && ...
                ~contains(n, 'q1') && ~contains(n, 'q2') && ~contains(n, 'q3') && ~contains(n, 'q4')

            score = 1;
            if contains(n, 'pupil'), score = score + 2; end

            if score > bestScore
                bestScore = score;
                bestVar = vars{i};
            end
        end
    end

    v = bestVar;

end

function v = find_var_exact_or_normalized(vars, aliases)

    v = '';

    % Exact, case-insensitive match
    for i = 1:numel(aliases)
        idx = find(strcmpi(vars, aliases{i}), 1, 'first');
        if ~isempty(idx)
            v = vars{idx};
            return;
        end
    end

    % Normalized match
    varsNorm = cellfun(@norm_name, vars, 'UniformOutput', false);

    for i = 1:numel(aliases)
        aNorm = norm_name(aliases{i});
        idx = find(strcmp(varsNorm, aNorm), 1, 'first');
        if ~isempty(idx)
            v = vars{idx};
            return;
        end
    end

end

function ST = run_corr_with_permutation(x, y, subj, corrType, doPermutationP, nPerm)

    x = as_numeric(x);
    y = as_numeric(y);
    subj = string(subj);

    ok = isfinite(x) & isfinite(y);

    x = x(ok);
    y = y(ok);
    subj = subj(ok);

    ST = struct();
    ST.x = x;
    ST.y = y;
    ST.subjects = subj;
    ST.N = numel(x);
    ST.r = NaN;
    ST.p = NaN;
    ST.pPerm = NaN;
    ST.CI95 = [NaN NaN];

    if ST.N < 5
        warning('Fewer than 5 complete cases. Correlation skipped.');
        return;
    end

    if numel(unique(x)) < 2 || numel(unique(y)) < 2
        warning('X or Y has fewer than 2 unique values. Correlation skipped.');
        return;
    end

    [r, p] = corr(x, y, 'Type', corrType, 'Rows', 'complete');

    ST.r = r;
    ST.p = p;

    if strcmpi(corrType, 'Pearson')
        ST.CI95 = pearson_ci(r, ST.N, 0.05);
    end

    if doPermutationP
        rNull = nan(nPerm,1);

        for ip = 1:nPerm
            yp = y(randperm(numel(y)));
            rNull(ip) = corr(x, yp, 'Type', corrType, 'Rows', 'complete');
        end

        ST.pPerm = (sum(abs(rNull) >= abs(r)) + 1) / (nPerm + 1);
    end

end

function plot_corr_axis(ax, x, y, subj, ST, titleStr, xLabelStr, yLabelStr, ...
                        corrType, showSubjectLabels, jitterY_forPlotOnly, jitterAmount)

    axes(ax); %#ok<LAXES>
    hold(ax, 'on');

    yPlot = y;

    if jitterY_forPlotOnly
        yPlot = y + (rand(size(y)) - 0.5) * 2 * jitterAmount;
    end

    scatter(ax, x, yPlot, 82, 'filled', ...
        'MarkerFaceColor', [0.20 0.20 0.20], ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.82);

    % Regression line is visualization only.
    % Statistics are Spearman/Pearson on the non-jittered true data.
    if numel(unique(x)) > 1 && numel(unique(y)) > 1
        xx = linspace(min(x), max(x), 100);
        b = polyfit(x, y, 1);
        yy = polyval(b, xx);
        plot(ax, xx, yy, 'k-', 'LineWidth', 2.4);
    end

    xlabel(ax, xLabelStr, 'FontSize', 13, 'FontWeight','bold');
    ylabel(ax, yLabelStr, 'FontSize', 13, 'FontWeight','bold');
    title(ax, titleStr, 'FontSize', 14, 'FontWeight','bold');

    box(ax, 'off');
    set(ax, 'FontSize', 12, 'LineWidth', 1.25, 'TickDir','out');
    grid(ax, 'on');

    xlim(ax, pad_limits(x, 0.12));

    if min(y) >= 0 && max(y) <= 2
        ylim(ax, [-0.25 2.25]);
        yticks(ax, 0:1:2);
    elseif min(y) >= 0 && max(y) <= 4
        ylim(ax, [-0.35 4.35]);
        yticks(ax, 0:1:4);
    else
        ylim(ax, pad_limits(y, 0.18));
    end

    if strcmpi(corrType, 'Pearson')
        statTxt = sprintf('N = %d\nr = %.2f\np = %.3f\nperm p = %.3f\np_{Holm} = %.3f', ...
            ST.N, ST.r, ST.p, ST.pPerm, ST.pHolm);
    else
        statTxt = sprintf('N = %d\n\\rho = %.2f\np = %.3f\nperm p = %.3f\np_{Holm} = %.3f', ...
            ST.N, ST.r, ST.p, ST.pPerm, ST.pHolm);
    end

    xl = xlim(ax);
    yl = ylim(ax);

    text(ax, xl(1) + 0.045*range(xl), yl(2) - 0.07*range(yl), statTxt, ...
        'FontSize', 11.5, ...
        'VerticalAlignment','top', ...
        'BackgroundColor','w', ...
        'EdgeColor',[0.80 0.80 0.80], ...
        'Margin',6);

    if showSubjectLabels
        for i = 1:numel(x)
            text(ax, x(i), yPlot(i), ['  ' char(subj(i))], ...
                'FontSize', 8, ...
                'Color', [0.25 0.25 0.25]);
        end
    end

end

function pAdj = holm_adjust(p)

    p = p(:);
    pAdj = nan(size(p));

    ok = isfinite(p);
    pOk = p(ok);

    if isempty(pOk)
        return;
    end

    [pSort, order] = sort(pOk);
    m = numel(pSort);

    adjSort = (m - (1:m)' + 1) .* pSort;
    adjSort = cummax(adjSort);
    adjSort(adjSort > 1) = 1;

    temp = nan(size(pOk));
    temp(order) = adjSort;

    pAdj(ok) = temp;

end

function ci = pearson_ci(r, n, alpha)

    if n <= 3 || abs(r) >= 1 || ~isfinite(r)
        ci = [NaN NaN];
        return;
    end

    z = atanh(r);
    se = 1 / sqrt(n - 3);

    % zcrit = norminv(1 - alpha/2), without requiring norminv
    zcrit = sqrt(2) * erfinv(1 - alpha);

    ci = tanh([z - zcrit*se, z + zcrit*se]);

end

function lim = pad_limits(x, frac)

    x = x(isfinite(x));

    if isempty(x)
        lim = [0 1];
        return;
    end

    if min(x) == max(x)
        lim = [min(x)-1, max(x)+1];
        return;
    end

    pad = frac * range(x);
    lim = [min(x)-pad, max(x)+pad];

end

function x = as_numeric(x)

    if isnumeric(x)
        x = double(x(:));
    elseif islogical(x)
        x = double(x(:));
    elseif isstring(x)
        x = str2double(x(:));
    elseif iscell(x)
        x = str2double(string(x(:)));
    elseif iscategorical(x)
        x = str2double(string(x(:)));
    else
        try
            x = double(x(:));
        catch
            x = str2double(string(x(:)));
        end
    end

end

function tf = any_contains(str, tokenCell)

    tf = false;

    for i = 1:numel(tokenCell)
        if contains(str, tokenCell{i})
            tf = true;
            return;
        end
    end

end

function n = norm_name(s)

    s = char(string(s));
    n = lower(regexprep(s, '[^a-zA-Z0-9]', ''));

end

function lbl = make_pretty_label(v)

    lbl = string(v);

    lbl = replace(lbl, '_', ' ');
    lbl = replace(lbl, 'E1', 'Counting (E1)');
    lbl = replace(lbl, 'E2', 'Visual interference (E2)');
    lbl = replace(lbl, 'RAS', 'RAS');
    lbl = replace(lbl, 'global effect', 'global pupil effect');
    lbl = replace(lbl, 'Global', 'global');
    lbl = replace(lbl, 'fromQ', '');
    lbl = replace(lbl, 'Q3', 'Q3');
    lbl = replace(lbl, 'Q4', 'Q4');

    lbl = char(lbl);

end