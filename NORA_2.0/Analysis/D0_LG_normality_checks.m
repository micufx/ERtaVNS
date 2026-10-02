clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% LG NORMALITY CHECKS: LILLIEFORS TESTS + Q-Q PLOTS
%
% Corrected standalone diagnostic script for the secondary parametric
% analyses in the local-global EEG-pupillometry study.
%
% IMPORTANT:
%   - This script reads participant-level summary tables already exported
%     by the pupil and EEG statistics pipelines.
%   - It does NOT repeat preprocessing.
%   - It does NOT modify the original statistical results.
%
% It evaluates normality of:
%   1) paired E1-E2 difference scores for Standard and Deviant conditions;
%   2) within-task Deviant-Standard difference scores;
%   3) the Task, Condition, and Task x Condition contrast scores underlying
%      each 2 x 2 repeated-measures ANOVA.
%
% Requires:
%   Statistics and Machine Learning Toolbox (lillietest, qqplot).
%
% Miguel Contreras-Altamirano, 2026

%% ========================= USER SETTINGS ================================

grandFolder = 'L:\Downloads\scripts_Vagus\matlab\testing\out\GRAND_AVG';

alphaNormality = 0.05;
makeQQPlots    = true;

% Set to false if slope analyses are not part of the manuscript.
includeSlope = true;

checkPupil = true;
checkEEG   = true;

%% ========================= BASIC CHECKS ================================

if ~exist(grandFolder, 'dir')
    error('GRAND_AVG folder not found: %s', grandFolder);
end

if exist('lillietest', 'file') ~= 2
    error(['lillietest was not found. This script requires the MATLAB ' ...
           'Statistics and Machine Learning Toolbox.']);
end

if exist('qqplot', 'file') ~= 2 && makeQQPlots
    warning('qqplot was not found. Q-Q plots will be disabled.');
    makeQQPlots = false;
end

outFolder = fullfile(grandFolder, 'NORMALITY_CHECKS_LILLIEFORS');
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

ALL_RESULTS = table();
filesUsed   = strings(0,1);

%% ========================= FIND SOURCE TABLES ===========================

% Recursively search GRAND_AVG and its subfolders. This is intentionally
% broader than the previous version so tables inside folders such as
% GRAND_AVG\tables are found automatically.

allCsv = dir(fullfile(grandFolder, '**', '*.csv'));

if isempty(allCsv)
    error('No CSV files were found under: %s', grandFolder);
end

allPaths = strings(numel(allCsv),1);
allNames = strings(numel(allCsv),1);

for i = 1:numel(allCsv)
    allPaths(i) = string(fullfile(allCsv(i).folder, allCsv(i).name));
    allNames(i) = string(allCsv(i).name);
end

%% ========================= PUPIL DATA =================================

if checkPupil

    pupilMask = contains(allNames, "PUPIL", 'IgnoreCase', true) & ...
                contains(allNames, "subject", 'IgnoreCase', true) & ...
                contains(allNames, "value", 'IgnoreCase', true);

    pupilCandidates = allPaths(pupilMask);

    % Prefer the ANOVA-specific table, then the general RainBox table.
    pupilFile = choose_best_table( ...
        pupilCandidates, ...
        ["PUPIL_RAINBOX_COMPARE_ANOVA_subject_values.csv", ...
         "PUPIL_RAINBOX_COMPARE_subject_values.csv"]);

    if strlength(pupilFile) == 0
        warning(['No pupil subject-values CSV was found recursively under GRAND_AVG. ' ...
                 'Expected a file such as PUPIL_RAINBOX_COMPARE_subject_values.csv.']);
    else
        fprintf('\n============================================================\n');
        fprintf('PUPIL SOURCE TABLE\n%s\n', pupilFile);
        fprintf('============================================================\n');

        T = readtable(pupilFile, 'VariableNamingRule', 'preserve');

        validate_subject_table(T, "Pupil", pupilFile);

        R = run_normality_family_checks( ...
            T, "Pupil", "Pupil", alphaNormality, includeSlope, ...
            makeQQPlots, outFolder);

        ALL_RESULTS = append_table(ALL_RESULTS, R);
        filesUsed(end+1,1) = pupilFile; %#ok<SAGROW>
    end
end

%% ========================= EEG DATA ===================================

if checkEEG

    eegMask = contains(allNames, "EEG", 'IgnoreCase', true) & ...
              contains(allNames, "subject", 'IgnoreCase', true) & ...
              contains(allNames, "value", 'IgnoreCase', true);

    eegCandidates = allPaths(eegMask);

    % Prefer direct RainBox exports because these are the participant-level
    % values used by the between-task EEG analyses.
    preferredEEG = eegCandidates( ...
        contains(eegCandidates, "EEG_RAINBOX_COMPARE_", 'IgnoreCase', true));

    if isempty(preferredEEG)
        preferredEEG = eegCandidates;
    end

    % Remove obvious output tables from a previous normality run.
    preferredEEG = preferredEEG( ...
        ~contains(preferredEEG, "Lilliefors", 'IgnoreCase', true) & ...
        ~contains(preferredEEG, "normality", 'IgnoreCase', true));

    if isempty(preferredEEG)
        warning(['No EEG subject-values CSV was found recursively under GRAND_AVG. ' ...
                 'Expected a file such as EEG_RAINBOX_COMPARE_Pz_subject_values.csv.']);
    else

        % Avoid processing duplicate representations of the same source.
        preferredEEG = unique(preferredEEG, 'stable');

        for iFile = 1:numel(preferredEEG)

            eegFile = preferredEEG(iFile);

            if any(filesUsed == eegFile)
                continue;
            end

            sourceLabel = infer_eeg_source_label(eegFile);

            fprintf('\n============================================================\n');
            fprintf('EEG SOURCE TABLE (%s)\n%s\n', sourceLabel, eegFile);
            fprintf('============================================================\n');

            T = readtable(eegFile, 'VariableNamingRule', 'preserve');

            % Skip unrelated EEG tables that happen to contain "subject"
            % and "values" but do not contain the required four cells.
            if ~has_required_subject_columns(T)
                fprintf('Skipping unrelated EEG table: %s\n', eegFile);
                continue;
            end

            R = run_normality_family_checks( ...
                T, "EEG", sourceLabel, alphaNormality, includeSlope, ...
                makeQQPlots, outFolder);

            ALL_RESULTS = append_table(ALL_RESULTS, R);
            filesUsed(end+1,1) = eegFile; %#ok<SAGROW>
        end
    end
end

%% ========================= SAVE OUTPUTS ================================

if isempty(ALL_RESULTS)
    error(['No normality analyses were completed. Check that the participant-level ' ...
           'subject-values CSV files exist and contain the expected columns.']);
end

outCsv  = fullfile(outFolder, 'LG_Lilliefors_normality_all_results.csv');
outXlsx = fullfile(outFolder, 'LG_Lilliefors_normality_all_results.xlsx');
outMat  = fullfile(outFolder, 'LG_Lilliefors_normality_all_results.mat');

writetable(ALL_RESULTS, outCsv);
writetable(ALL_RESULTS, outXlsx, 'FileType', 'spreadsheet');
save(outMat, 'ALL_RESULTS', 'filesUsed', 'alphaNormality');

FLAGGED = ALL_RESULTS(ALL_RESULTS.RejectNormality == 1, :);
outFlagged = fullfile(outFolder, 'LG_Lilliefors_normality_flagged_only.csv');
writetable(FLAGGED, outFlagged);

% Modality summary
SUMMARY = groupsummary( ...
    ALL_RESULTS, ...
    {'Modality','Source'}, ...
    {'sum'}, ...
    'RejectNormality');

SUMMARY.Properties.VariableNames( ...
    strcmp(SUMMARY.Properties.VariableNames,'sum_RejectNormality')) = ...
    {'NumberRejected'};

outSummary = fullfile(outFolder, 'LG_Lilliefors_normality_summary_by_source.csv');
writetable(SUMMARY, outSummary);

fprintf('\n============================================================\n');
fprintf('LILLIEFORS NORMALITY CHECKS COMPLETED\n');
fprintf('Alpha = %.3f\n', alphaNormality);
fprintf('Source files used: %d\n', numel(filesUsed));

for i = 1:numel(filesUsed)
    fprintf('  %s\n', filesUsed(i));
end

fprintf('\nTestable distributions: %d\n', sum(isfinite(ALL_RESULTS.PValue)));
fprintf('Normality not rejected: %d\n', sum(ALL_RESULTS.RejectNormality == 0));
fprintf('Normality rejected:     %d\n', sum(ALL_RESULTS.RejectNormality == 1));
fprintf('Not testable:           %d\n', sum(~isfinite(ALL_RESULTS.PValue)));

fprintf('\nSaved results:\n');
fprintf('%s\n%s\n%s\n%s\n', outCsv, outXlsx, outFlagged, outSummary);
fprintf('Q-Q plots folder:\n%s\n', outFolder);
fprintf('============================================================\n\n');

disp(ALL_RESULTS(:, { ...
    'Modality','Source','Effect','Measure','Distribution', ...
    'N','PValue','RejectNormality','Interpretation'}));

%% ========================= LOCAL FUNCTIONS =============================

function T = append_table(T, newRows)

if isempty(newRows)
    return;
end

if isempty(T)
    T = newRows;
else
    T = [T; newRows];
end
end

function chosen = choose_best_table(candidates, preferredNames)

chosen = "";

if isempty(candidates)
    return;
end

% First try exact preferred filenames in the requested order.
for iPref = 1:numel(preferredNames)
    idx = find(endsWith(candidates, preferredNames(iPref), ...
        'IgnoreCase', true), 1);

    if ~isempty(idx)
        chosen = candidates(idx);
        return;
    end
end

% Otherwise use the first recursively discovered candidate.
chosen = candidates(1);
end

function sourceLabel = infer_eeg_source_label(filePath)

s = string(filePath);

if contains(s, "Pz", 'IgnoreCase', true)
    sourceLabel = "EEG Pz";
elseif contains(s, "Fz", 'IgnoreCase', true)
    sourceLabel = "EEG Fz";
elseif contains(s, "Cz", 'IgnoreCase', true)
    sourceLabel = "EEG Cz";
elseif contains(s, "ROI", 'IgnoreCase', true)
    sourceLabel = "EEG ROI";
else
    sourceLabel = "EEG";
end
end

function validate_subject_table(T, sourceName, filePath)

if ~has_required_subject_columns(T)
    error(['The selected %s table does not contain the expected Local/Global, ' ...
           'Standard/Deviant, E1/E2 participant-level columns:\n%s'], ...
           sourceName, filePath);
end
end

function tf = has_required_subject_columns(T)

vars = string(T.Properties.VariableNames);

requiredGroups = { ...
    ["LocalStandard_Mean_E1","LocalStandard_Mean_E1_uV"], ...
    ["LocalStandard_Mean_E2","LocalStandard_Mean_E2_uV"], ...
    ["LocalDeviant_Mean_E1","LocalDeviant_Mean_E1_uV"], ...
    ["LocalDeviant_Mean_E2","LocalDeviant_Mean_E2_uV"], ...
    ["GlobalStandard_Mean_E1","GlobalStandard_Mean_E1_uV"], ...
    ["GlobalStandard_Mean_E2","GlobalStandard_Mean_E2_uV"], ...
    ["GlobalDeviant_Mean_E1","GlobalDeviant_Mean_E1_uV"], ...
    ["GlobalDeviant_Mean_E2","GlobalDeviant_Mean_E2_uV"]};

tf = true;

for i = 1:numel(requiredGroups)
    found = any(arrayfun(@(x) any(strcmpi(vars, x)), requiredGroups{i}));
    if ~found
        tf = false;
        return;
    end
end
end

function RESULTS = run_normality_family_checks( ...
    T, modality, sourceLabel, alphaNormality, includeSlope, ...
    makeQQPlots, outFolder)

RESULTS = table();

effects = ["Local","Global"];

measureNames = ["Peak latency","Mean response"];
measureKeys  = ["Latency","Mean"];

if includeSlope
    measureNames(end+1) = "Slope";
    measureKeys(end+1)  = "Slope";
end

for iEffect = 1:numel(effects)

    effectName = effects(iEffect);

    for iMeasure = 1:numel(measureKeys)

        measureName = measureNames(iMeasure);
        measureKey  = measureKeys(iMeasure);

        [stdE1Name, stdE2Name, devE1Name, devE2Name] = ...
            resolve_four_columns(T, effectName, measureKey);

        if strlength(stdE1Name) == 0
            fprintf('Skipping unavailable family: %s | %s | %s\n', ...
                sourceLabel, effectName, measureName);
            continue;
        end

        stdE1 = double(T.(char(stdE1Name)));
        stdE2 = double(T.(char(stdE2Name)));
        devE1 = double(T.(char(devE1Name)));
        devE2 = double(T.(char(devE2Name)));

        distNames = [ ...
            "ANOVA Task contrast", ...
            "ANOVA Condition contrast", ...
            "ANOVA Task x Condition contrast", ...
            "Paired E1-E2 difference | Standard", ...
            "Paired E1-E2 difference | Deviant", ...
            "Within E1 difference | Deviant-Standard", ...
            "Within E2 difference | Deviant-Standard"];

        distValues = cell(1, numel(distNames));

        complete = isfinite(stdE1) & isfinite(devE1) & ...
                   isfinite(stdE2) & isfinite(devE2);

        yS1 = stdE1(complete);
        yD1 = devE1(complete);
        yS2 = stdE2(complete);
        yD2 = devE2(complete);

        % Exact subject-level contrasts underlying a 2 x 2 within-subject
        % ANOVA with factors Task and Condition.
        distValues{1} = ((yS1 + yD1) ./ 2) - ((yS2 + yD2) ./ 2);
        distValues{2} = ((yD1 + yD2) ./ 2) - ((yS1 + yS2) ./ 2);
        distValues{3} = (yD1 - yS1) - (yD2 - yS2);

        % Paired and within-task difference scores.
        distValues{4} = paired_difference(stdE1, stdE2);
        distValues{5} = paired_difference(devE1, devE2);
        distValues{6} = paired_difference(devE1, stdE1);
        distValues{7} = paired_difference(devE2, stdE2);

        familyRows = table();

        if makeQQPlots
            fig = figure( ...
                'Color','w', ...
                'Units','normalized', ...
                'OuterPosition',[0 0 1 1], ...
                'Visible','off');

            tiledlayout(fig, 2, 4, ...
                'TileSpacing','compact', ...
                'Padding','compact');
        else
            fig = [];
        end

        for iDist = 1:numel(distNames)

            x = double(distValues{iDist}(:));
            x = x(isfinite(x));

            [row, plotInfo] = run_one_lilliefors( ...
                x, modality, sourceLabel, effectName, measureName, ...
                distNames(iDist), alphaNormality);

            familyRows = append_table(familyRows, row);

            if makeQQPlots
                ax = nexttile;
                plot_qq_panel(ax, x, distNames(iDist), plotInfo);
            end
        end

        if makeQQPlots

            ax = nexttile;
            axis(ax, 'off');

            text(ax, 0, 0.90, sprintf('%s | %s | %s', ...
                sourceLabel, effectName, measureName), ...
                'FontWeight','bold', 'Interpreter','none');

            text(ax, 0, 0.65, sprintf('Lilliefors alpha = %.3f', ...
                alphaNormality), 'Interpreter','none');

            text(ax, 0, 0.40, ...
                'h = 0: normality not rejected', ...
                'Interpreter','none');

            text(ax, 0, 0.20, ...
                'h = 1: normality rejected', ...
                'Interpreter','none');

            sgtitle(fig, sprintf( ...
                'Normality diagnostics | %s | %s | %s', ...
                sourceLabel, effectName, measureName), ...
                'FontWeight','bold', 'Interpreter','none');

            safeName = matlab.lang.makeValidName( ...
                char(sourceLabel + "_" + effectName + "_" + measureName));

            outPng = fullfile(outFolder, ['QQ_' safeName '.png']);

            exportgraphics(fig, outPng, ...
                'Resolution', 200, ...
                'BackgroundColor','white');

            close(fig);
        end

        RESULTS = append_table(RESULTS, familyRows);
    end
end
end

function d = paired_difference(a, b)

keep = isfinite(a) & isfinite(b);
d = double(a(keep) - b(keep));
d = d(:);
end

function [stdE1Name, stdE2Name, devE1Name, devE2Name] = ...
    resolve_four_columns(T, effectName, measureKey)

vars = string(T.Properties.VariableNames);

stdBase = effectName + "Standard";
devBase = effectName + "Deviant";

switch measureKey

    case "Latency"
        stdE1Candidates = stdBase + ["_Latency_E1_ms"];
        stdE2Candidates = stdBase + ["_Latency_E2_ms"];
        devE1Candidates = devBase + ["_Latency_E1_ms"];
        devE2Candidates = devBase + ["_Latency_E2_ms"];

    case "Mean"
        % Pupil tables use no unit suffix; EEG tables use _uV.
        stdE1Candidates = stdBase + ["_Mean_E1", "_Mean_E1_uV"];
        stdE2Candidates = stdBase + ["_Mean_E2", "_Mean_E2_uV"];
        devE1Candidates = devBase + ["_Mean_E1", "_Mean_E1_uV"];
        devE2Candidates = devBase + ["_Mean_E2", "_Mean_E2_uV"];

    case "Slope"
        stdE1Candidates = stdBase + ["_Slope_E1"];
        stdE2Candidates = stdBase + ["_Slope_E2"];
        devE1Candidates = devBase + ["_Slope_E1"];
        devE2Candidates = devBase + ["_Slope_E2"];

    otherwise
        error('Unknown measure key: %s', measureKey);
end

stdE1Name = find_column(vars, stdE1Candidates);
stdE2Name = find_column(vars, stdE2Candidates);
devE1Name = find_column(vars, devE1Candidates);
devE2Name = find_column(vars, devE2Candidates);

if any(strlength([stdE1Name,stdE2Name,devE1Name,devE2Name]) == 0)
    stdE1Name = "";
    stdE2Name = "";
    devE1Name = "";
    devE2Name = "";
end
end

function found = find_column(vars, candidates)

found = "";

for i = 1:numel(candidates)
    idx = find(strcmpi(vars, candidates(i)), 1);

    if ~isempty(idx)
        found = vars(idx);
        return;
    end
end
end

function [row, info] = run_one_lilliefors( ...
    x, modality, sourceLabel, effectName, measureName, ...
    distributionName, alphaNormality)

n = numel(x);

h       = NaN;
p       = NaN;
kstat   = NaN;
critVal = NaN;

if n < 4
    interpretation = "Not testable: fewer than four valid observations";

elseif std(x, 0) <= eps(max(abs(x)) + 1)
    interpretation = "Not testable: distribution is constant";

else
    try
        [h, p, kstat, critVal] = lillietest( ...
            x, ...
            'Alpha', alphaNormality, ...
            'Distribution', 'normal');

        if h == 0
            interpretation = "Normality not rejected";
        else
            interpretation = "Normality rejected";
        end

    catch
        interpretation = "Not testable: lillietest returned an error";
        h = NaN;
        p = NaN;
        kstat = NaN;
        critVal = NaN;
    end
end

row = table( ...
    string(modality), ...
    string(sourceLabel), ...
    string(effectName), ...
    string(measureName), ...
    string(distributionName), ...
    n, ...
    mean(x, 'omitnan'), ...
    std(x, 0, 'omitnan'), ...
    skewness_safe(x), ...
    kurtosis_safe(x), ...
    kstat, ...
    critVal, ...
    p, ...
    h, ...
    string(interpretation), ...
    'VariableNames', { ...
    'Modality', ...
    'Source', ...
    'Effect', ...
    'Measure', ...
    'Distribution', ...
    'N', ...
    'Mean', ...
    'SD', ...
    'Skewness', ...
    'Kurtosis', ...
    'LillieforsStatistic', ...
    'CriticalValue', ...
    'PValue', ...
    'RejectNormality', ...
    'Interpretation'});

info = struct();
info.h = h;
info.p = p;
info.n = n;
info.interpretation = interpretation;
end

function plot_qq_panel(ax, x, panelTitle, info)

if numel(x) < 4 || std(x, 0) <= eps(max(abs(x)) + 1)
    axis(ax, 'off');

    text(ax, 0.05, 0.60, panelTitle, ...
        'FontWeight','bold', ...
        'Interpreter','none');

    text(ax, 0.05, 0.35, info.interpretation, ...
        'Interpreter','none');

    return;
end

axes(ax); %#ok<LAXES>
qqplot(x);
grid(ax, 'on');

if isfinite(info.p)
    title(ax, sprintf('%s\np = %.4g | h = %d | n = %d', ...
        panelTitle, info.p, info.h, info.n), ...
        'Interpreter','none', ...
        'FontSize', 9);
else
    title(ax, sprintf('%s\nnot testable | n = %d', ...
        panelTitle, info.n), ...
        'Interpreter','none', ...
        'FontSize', 9);
end
end

function s = skewness_safe(x)

if numel(x) < 3 || std(x,0) == 0
    s = NaN;
else
    s = skewness(x, 0);
end
end

function k = kurtosis_safe(x)

if numel(x) < 4 || std(x,0) == 0
    k = NaN;
else
    k = kurtosis(x, 0);
end
end
