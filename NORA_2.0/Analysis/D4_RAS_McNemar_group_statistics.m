clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% D3_RAS_McNemar_group_statistics.m
%
% PURPOSE
%   Stand-alone group-level analysis of behavioral Regularity Awareness
%   Score (RAS) data across the two within-subject tasks:
%
%       E1 = Counting
%       E2 = Visual interference
%
% PRIMARY ANALYSIS
%   Total RAS is computed from Q1-Q4 for each experiment (range 0-4).
%   Participants are classified as:
%
%       Maximum awareness       = RAS 4/4
%       Below-maximum awareness = RAS < 4/4
%
%   Because the SAME participants completed E1 and E2, the primary test is
%   a two-sided EXACT McNemar test based on the discordant paired outcomes.
%
%   No multiple-comparison correction is applied to this single,
%   prespecified primary behavioral comparison.
%
% SENSITIVITY ANALYSIS
%   For transparency, the script also calculates a two-sided Fisher exact
%   test from the marginal 2x2 task-by-awareness table. Fisher's test treats
%   the two task samples as independent and is therefore NOT the preferred
%   inferential test for this within-subject design.
%
% OPTIONAL EXPLORATORY ANALYSES
%   Q3 and Q4 are compared separately between E1 and E2 using exact
%   McNemar tests. Holm correction is applied across these two item-level
%   tests. These analyses are exploratory and need not be reported.
%
% INPUT
%   The script uses the same Pupil_RAS.mat file and Combined table/struct
%   used by D1_Pupil_RAS_correlation.m.
%
% OUTPUT
%   - Participant-level diagnostic CSV
%   - 2x2 paired transition table CSV
%   - Primary and sensitivity results CSV
%   - Optional item-level results CSV
%   - MATLAB results MAT file
%   - Plain-text report
%   - Publication-style paired transition figure:
%       PNG (300 dpi), PDF (vector), and FIG
%
% Miguel Contreras-Altamirano, 2026

%% ===================== USER SETTINGS =====================

cfg = struct();

% Same input location used in D1_Pupil_RAS_correlation.m
cfg.basepath = 'L:\Downloads\scripts_Vagus\matlab\testing\out\GRAND_AVG\tables\Pupil_RAS';
cfg.matFile  = fullfile(cfg.basepath, 'Pupil_RAS.mat');

% Independent output folder for this behavioral analysis
cfg.outFolder = fullfile(cfg.basepath, 'RAS_McNemar_group_statistics');

% % Apply the same exclusions used in the physiological analyses.
% % If these subjects are already absent from Pupil_RAS.mat, nothing changes.
% cfg.excludeSubjects = ["S21","S12","S23","S02","S29"];
cfg.excludeSubjects = strings(0,1);

% Optional exploratory item-level tests (Q3 and Q4)
cfg.runExploratoryItemTests = true;

% Figure settings
cfg.saveFigure = true;
cfg.exportDPI  = 300;

if ~exist(cfg.outFolder, 'dir')
    mkdir(cfg.outFolder);
end

%% ===================== LOAD DATA =====================

if exist(cfg.matFile, 'file') ~= 2
    error('Input MAT file not found:\n%s', cfg.matFile);
end

S = load(cfg.matFile);
D = load_combined_as_table(S);
D = ensure_subject_column(D);
D = add_ras_scores(D);

fprintf('\nLoaded %d rows from:\n  %s\n', height(D), cfg.matFile);

%% ===================== APPLY EXCLUSIONS =====================

if ~isempty(cfg.excludeSubjects)
    beforeN = height(D);
    D = D(~ismember(string(D.Subject), cfg.excludeSubjects), :);
    fprintf('Excluded %d row(s) using the configured exclusion list.\n', ...
        beforeN - height(D));
end

% One row per participant is required.
[uSub, ~, iSub] = unique(string(D.Subject));
subCounts = accumarray(iSub, 1);
if any(subCounts > 1)
    duplicated = uSub(subCounts > 1);
    error('Duplicate participant rows found: %s', strjoin(duplicated, ', '));
end

%% ===================== CHECK AND PREPARE RAS DATA =====================

requiredVars = { ...
    'Q1_E1','Q2_E1','Q3_E1','Q4_E1', ...
    'Q1_E2','Q2_E2','Q3_E2','Q4_E2'};

assert_variables_exist(D, requiredVars);

% Convert all questionnaire variables to numeric.
for i = 1:numel(requiredVars)
    D.(requiredVars{i}) = as_numeric(D.(requiredVars{i}));
end

% Keep complete paired questionnaire records.
Qall = D{:, requiredVars};
completeMask = all(isfinite(Qall), 2);

if any(~completeMask)
    fprintf('Removed %d participant(s) with incomplete E1/E2 RAS data.\n', ...
        sum(~completeMask));
end

D = D(completeMask, :);

% Verify binary item coding.
Qall = D{:, requiredVars};
badValueMask = ~(Qall == 0 | Qall == 1);
if any(badValueMask, 'all')
    badValues = unique(Qall(badValueMask));
    error(['RAS item values must be coded 0/1. Unexpected value(s): %s'], ...
        strjoin(string(badValues(:)'), ', '));
end

% Total and global scores.
D.RAS_E1_Total  = D.Q1_E1 + D.Q2_E1 + D.Q3_E1 + D.Q4_E1;
D.RAS_E2_Total  = D.Q1_E2 + D.Q2_E2 + D.Q3_E2 + D.Q4_E2;
D.RAS_E1_Global = D.Q3_E1 + D.Q4_E1;
D.RAS_E2_Global = D.Q3_E2 + D.Q4_E2;

% Primary binary endpoint.
D.E1_MaximumAwareness = D.RAS_E1_Total == 4;
D.E2_MaximumAwareness = D.RAS_E2_Total == 4;

N = height(D);
if N < 2
    error('Fewer than two complete paired participants remain.');
end

%% ===================== PAIRED TRANSITION COUNTS =====================

e1Max = logical(D.E1_MaximumAwareness);
e2Max = logical(D.E2_MaximumAwareness);

% Paired 2x2 table:
%                    E2 maximum   E2 below maximum
% E1 maximum            n11              n10
% E1 below maximum      n01              n00

n11 = sum( e1Max &  e2Max);  % remained maximum
n10 = sum( e1Max & ~e2Max);  % decreased
n01 = sum(~e1Max &  e2Max);  % increased
n00 = sum(~e1Max & ~e2Max);  % remained below maximum

nE1Max = sum(e1Max);
nE2Max = sum(e2Max);

pctE1Max = 100 * nE1Max / N;
pctE2Max = 100 * nE2Max / N;
differencePP = pctE2Max - pctE1Max; % E2 minus E1, percentage points

D.Transition = strings(N,1);
D.Transition( e1Max &  e2Max) = "Remained maximum";
D.Transition( e1Max & ~e2Max) = "Decreased";
D.Transition(~e1Max &  e2Max) = "Increased";
D.Transition(~e1Max & ~e2Max) = "Remained below maximum";

%% ===================== PRIMARY EXACT McNEMAR TEST =====================

% McNemar uses only discordant pairs:
%   b = maximum in E1, below maximum in E2
%   c = below maximum in E1, maximum in E2

b = n10;
c = n01;

pMcNemar = exact_mcnemar_two_sided(b, c);

if c == 0 && b > 0
    matchedPairsOR = Inf;
elseif b == 0 && c > 0
    matchedPairsOR = 0;
elseif b == 0 && c == 0
    matchedPairsOR = NaN;
else
    matchedPairsOR = b / c;
end

%% ===================== FISHER SENSITIVITY TEST =====================

% Marginal task-by-awareness table:
%                       Maximum   Below maximum
% Counting (E1)          nE1Max      N-nE1Max
% Visual interference    nE2Max      N-nE2Max
%
% This ignores pairing and is saved only because Fisher's test was
% specifically suggested. The exact McNemar test remains the primary test.

fisherTable = [ ...
    nE1Max, N - nE1Max; ...
    nE2Max, N - nE2Max];

[pFisher, fisherOR] = fisher_exact_two_sided(fisherTable);

%% ===================== SAVE SUBJECT-LEVEL DATA =====================

subjectVars = { ...
    'Subject', ...
    'Q1_E1','Q2_E1','Q3_E1','Q4_E1','RAS_E1_Global','RAS_E1_Total', ...
    'Q1_E2','Q2_E2','Q3_E2','Q4_E2','RAS_E2_Global','RAS_E2_Total', ...
    'E1_MaximumAwareness','E2_MaximumAwareness','Transition'};

Tsubject = D(:, subjectVars);

subjectCsv = fullfile(cfg.outFolder, 'RAS_McNemar_subject_data.csv');
writetable(Tsubject, subjectCsv);

%% ===================== SAVE PAIRED TRANSITION TABLE =====================

Ttransition = array2table( ...
    [n11 n10; n01 n00], ...
    'VariableNames', {'E2_Maximum','E2_BelowMaximum'}, ...
    'RowNames', {'E1_Maximum','E1_BelowMaximum'});

transitionCsv = fullfile(cfg.outFolder, 'RAS_McNemar_paired_transition_table.csv');
writetable(Ttransition, transitionCsv, 'WriteRowNames', true);

%% ===================== SAVE PRIMARY RESULTS =====================

Tresults = table( ...
    N, ...
    nE1Max, pctE1Max, ...
    nE2Max, pctE2Max, ...
    differencePP, ...
    n11, n10, n01, n00, ...
    b, c, ...
    pMcNemar, matchedPairsOR, ...
    pFisher, fisherOR, ...
    false, ...
    "Single prespecified primary behavioral comparison", ...
    'VariableNames', { ...
    'N', ...
    'E1_Maximum_N','E1_Maximum_Percent', ...
    'E2_Maximum_N','E2_Maximum_Percent', ...
    'E2_minus_E1_PercentagePoints', ...
    'RemainedMaximum_N','Decreased_N','Increased_N','RemainedBelow_N', ...
    'Discordant_E1Max_to_E2Below_N','Discordant_E1Below_to_E2Max_N', ...
    'ExactMcNemar_P','MatchedPairsOddsRatio', ...
    'FisherExactSensitivity_P','FisherOddsRatio', ...
    'PrimaryMultipleCorrectionApplied','PrimaryCorrectionReason'});

resultsCsv = fullfile(cfg.outFolder, 'RAS_McNemar_primary_results.csv');
writetable(Tresults, resultsCsv);

%% ===================== OPTIONAL ITEM-LEVEL TESTS =====================

Titems = table();

if cfg.runExploratoryItemTests

    itemNames = ["Q3","Q4"];
    rawP = nan(numel(itemNames),1);

    for i = 1:numel(itemNames)

        item = itemNames(i);
        e1 = D.(char(item + "_E1")) == 1;
        e2 = D.(char(item + "_E2")) == 1;

        stayedCorrect = sum( e1 &  e2);
        lostCorrect   = sum( e1 & ~e2);
        gainedCorrect = sum(~e1 &  e2);
        stayedIncorrect = sum(~e1 & ~e2);

        rawP(i) = exact_mcnemar_two_sided(lostCorrect, gainedCorrect);

        row = table( ...
            item, ...
            sum(e1), 100*sum(e1)/N, ...
            sum(e2), 100*sum(e2)/N, ...
            stayedCorrect, lostCorrect, gainedCorrect, stayedIncorrect, ...
            rawP(i), ...
            'VariableNames', { ...
            'Item','E1_Correct_N','E1_Correct_Percent', ...
            'E2_Correct_N','E2_Correct_Percent', ...
            'RemainedCorrect_N','LostCorrect_N','GainedCorrect_N', ...
            'RemainedIncorrect_N','ExactMcNemar_P'});

        Titems = [Titems; row]; %#ok<AGROW>
    end

    Titems.HolmAdjusted_P = holm_adjust(rawP);
    Titems.AnalysisStatus = repmat( ...
        "Exploratory item-level analysis; Holm corrected across Q3 and Q4", ...
        height(Titems), 1);

    itemCsv = fullfile(cfg.outFolder, 'RAS_McNemar_exploratory_Q3_Q4_results.csv');
    writetable(Titems, itemCsv);
end

%% ===================== PRINT RESULTS =====================

fprintf('\n============================================================\n');
fprintf('RAS GROUP-LEVEL BEHAVIORAL ANALYSIS\n');
fprintf('============================================================\n');
fprintf('Complete paired sample: N = %d\n\n', N);

fprintf('Maximum awareness (total RAS = 4/4):\n');
fprintf('  Counting (E1):            %d/%d (%.1f%%)\n', nE1Max, N, pctE1Max);
fprintf('  Visual interference (E2): %d/%d (%.1f%%)\n', nE2Max, N, pctE2Max);
fprintf('  E2 - E1 difference:       %.1f percentage points\n\n', differencePP);

fprintf('Participant transitions:\n');
fprintf('  Remained maximum:         %d\n', n11);
fprintf('  Decreased:                %d\n', n10);
fprintf('  Increased:                %d\n', n01);
fprintf('  Remained below maximum:   %d\n\n', n00);

fprintf('PRIMARY PAIRED TEST:\n');
fprintf('  Two-sided exact McNemar p = %.6f\n', pMcNemar);
fprintf('  Discordant pairs: b = %d, c = %d\n', b, c);
fprintf('  Matched-pairs odds ratio (b/c) = %s\n', format_number(matchedPairsOR));
fprintf('  Multiple-comparison correction: not applied\n');
fprintf('  Reason: one prespecified primary behavioral comparison.\n\n');

fprintf('FISHER SENSITIVITY CHECK (ignores within-subject pairing):\n');
fprintf('  Two-sided Fisher exact p = %.6f\n', pFisher);
fprintf('  Fisher odds ratio = %s\n\n', format_number(fisherOR));

if cfg.runExploratoryItemTests
    fprintf('EXPLORATORY ITEM-LEVEL TESTS (Holm corrected across Q3 and Q4):\n');
    for i = 1:height(Titems)
        fprintf('  %s: raw p = %.6f, Holm p = %.6f\n', ...
            Titems.Item(i), Titems.ExactMcNemar_P(i), ...
            Titems.HolmAdjusted_P(i));
    end
    fprintf('\n');
end

%% ===================== CREATE TRANSITION FIGURE =====================

fig = make_transition_figure( ...
    N, n11, n10, n01, n00, ...
    nE1Max, nE2Max, pctE1Max, pctE2Max, ...
    differencePP, pMcNemar);

if cfg.saveFigure

    pngFile = fullfile(cfg.outFolder, 'RAS_awareness_transition_E1_vs_E2.png');
    pdfFile = fullfile(cfg.outFolder, 'RAS_awareness_transition_E1_vs_E2.pdf');
    figFile = fullfile(cfg.outFolder, 'RAS_awareness_transition_E1_vs_E2.fig');

    exportgraphics(fig, pngFile, 'Resolution', cfg.exportDPI, ...
        'BackgroundColor', 'white');
    exportgraphics(fig, pdfFile, 'ContentType', 'vector', ...
        'BackgroundColor', 'white');
    savefig(fig, figFile);

    fprintf('Saved figure:\n');
    fprintf('  %s\n', pngFile);
    fprintf('  %s\n', pdfFile);
    fprintf('  %s\n\n', figFile);
end

%% ===================== SAVE MAT AND TEXT REPORT =====================

resultsMat = fullfile(cfg.outFolder, 'RAS_McNemar_group_statistics.mat');
save(resultsMat, 'cfg', 'D', 'Tsubject', 'Ttransition', ...
    'Tresults', 'Titems', 'fisherTable', '-v7.3');

reportFile = fullfile(cfg.outFolder, 'RAS_McNemar_report.txt');
fid = fopen(reportFile, 'w');

if fid == -1
    warning('Could not create report file: %s', reportFile);
else
    fprintf(fid, 'RAS GROUP-LEVEL BEHAVIORAL ANALYSIS\n');
    fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));

    fprintf(fid, 'Complete paired sample: N = %d\n\n', N);
    fprintf(fid, 'Maximum awareness was defined as total RAS = 4/4.\n');
    fprintf(fid, 'Counting (E1): %d/%d (%.1f%%)\n', nE1Max, N, pctE1Max);
    fprintf(fid, 'Visual interference (E2): %d/%d (%.1f%%)\n', ...
        nE2Max, N, pctE2Max);
    fprintf(fid, 'E2 - E1 difference: %.1f percentage points\n\n', ...
        differencePP);

    fprintf(fid, 'Transitions:\n');
    fprintf(fid, 'Remained maximum: %d\n', n11);
    fprintf(fid, 'Decreased: %d\n', n10);
    fprintf(fid, 'Increased: %d\n', n01);
    fprintf(fid, 'Remained below maximum: %d\n\n', n00);

    fprintf(fid, 'Primary two-sided exact McNemar test:\n');
    fprintf(fid, 'b = %d, c = %d, p = %.6f\n', b, c, pMcNemar);
    fprintf(fid, 'Matched-pairs odds ratio = %s\n', ...
        format_number(matchedPairsOR));
    fprintf(fid, ['No multiple-comparison correction was applied because ' ...
        'this was one prespecified primary behavioral comparison.\n\n']);

    fprintf(fid, 'Fisher exact sensitivity check (pairing ignored):\n');
    fprintf(fid, 'p = %.6f, odds ratio = %s\n\n', ...
        pFisher, format_number(fisherOR));

    if cfg.runExploratoryItemTests
        fprintf(fid, 'Exploratory item-level exact McNemar tests:\n');
        for i = 1:height(Titems)
            fprintf(fid, '%s: raw p = %.6f, Holm p = %.6f\n', ...
                Titems.Item(i), Titems.ExactMcNemar_P(i), ...
                Titems.HolmAdjusted_P(i));
        end
    end

    fclose(fid);
end

fprintf('Saved result files to:\n  %s\n', cfg.outFolder);
fprintf('\nDone.\n');

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function D = load_combined_as_table(S)

    if isfield(S, 'Combined')
        A = S.Combined;
        fprintf('Using variable S.Combined.\n');
    else
        names = fieldnames(S);
        A = [];
        selectedName = '';

        for i = 1:numel(names)
            candidate = S.(names{i});
            if istable(candidate) || isstruct(candidate)
                A = candidate;
                selectedName = names{i};
                break;
            end
        end

        if isempty(A)
            error('No table or struct was found inside the MAT file.');
        end

        fprintf('No Combined variable found. Using S.%s.\n', selectedName);
    end

    if istable(A)
        D = A;
    elseif isstruct(A)
        if numel(A) > 1
            D = struct2table(A);
        else
            D = scalar_struct_to_table(A);
        end
    else
        error('Loaded object is neither a table nor a struct.');
    end
end

function T = scalar_struct_to_table(A)

    fields = fieldnames(A);
    T = table();

    for i = 1:numel(fields)

        name = fields{i};
        value = A.(name);

        if ischar(value)
            value = string(cellstr(value));
        elseif iscell(value)
            try
                value = string(value);
            catch
                % Retain the original cell data if conversion is unsuitable.
            end
        elseif iscategorical(value)
            value = string(value);
        end

        if isrow(value) && ~isscalar(value)
            value = value(:);
        end

        T.(name) = value;
    end
end

function D = ensure_subject_column(D)

    vars = string(D.Properties.VariableNames);

    exactIdx = find(strcmpi(vars, "Subject"), 1, 'first');

    if ~isempty(exactIdx)
        subjectData = D.(char(vars(exactIdx)));
    else
        idx = find(contains(lower(vars), "subject") | ...
                   contains(lower(vars), "participant"), 1, 'first');

        if isempty(idx)
            warning(['No participant column was found. Subject labels are ' ...
                'being generated from row order.']);
            subjectData = (1:height(D))';
        else
            subjectData = D.(char(vars(idx)));
            fprintf('Using %s as the participant column.\n', vars(idx));
        end
    end

    D.Subject = standardize_subject_labels(subjectData);
end

function subject = standardize_subject_labels(input)

    if isnumeric(input)
        subject = strings(numel(input),1);
        for i = 1:numel(input)
            subject(i) = sprintf('S%02d', input(i));
        end
        return;
    end

    input = string(input);
    subject = strings(numel(input),1);

    for i = 1:numel(input)

        current = strtrim(input(i));
        token = regexp(char(current), '\d+', 'match');

        if isempty(token)
            subject(i) = current;
        else
            number = str2double(token{end});
            if isfinite(number)
                subject(i) = sprintf('S%02d', number);
            else
                subject(i) = current;
            end
        end
    end
end

function D = add_ras_scores(D)

    questions = ["Q1","Q2","Q3","Q4"];
    experiments = ["E1","E2"];

    for e = 1:numel(experiments)
        for q = 1:numel(questions)

            canonical = questions(q) + "_" + experiments(e);
            sourceName = find_ras_question_variable( ...
                D, experiments(e), questions(q));

            if strlength(sourceName) == 0
                fprintf('\nAvailable variables:\n');
                disp(D.Properties.VariableNames');
                error('Could not find the RAS variable for %s %s.', ...
                    experiments(e), questions(q));
            end

            D.(char(canonical)) = as_numeric(D.(char(sourceName)));

            if sourceName ~= canonical
                fprintf('Mapped %-10s <- %s\n', canonical, sourceName);
            end
        end
    end
end

function variableName = find_ras_question_variable(D, experiment, question)

    vars = string(D.Properties.VariableNames);
    normalizedVars = normalize_name(vars);

    aliases = [ ...
        question + "_" + experiment, ...
        experiment + "_" + question, ...
        question + experiment, ...
        experiment + question];

    normalizedAliases = normalize_name(aliases);

    variableName = "";

    for i = 1:numel(normalizedAliases)
        idx = find(normalizedVars == normalizedAliases(i), 1, 'first');
        if ~isempty(idx)
            variableName = vars(idx);
            return;
        end
    end

    % Fallback: require both experiment and question tokens.
    expToken = lower(experiment);
    qToken = lower(question);

    candidates = contains(normalizedVars, expToken) & ...
                 contains(normalizedVars, qToken);

    idx = find(candidates, 1, 'first');
    if ~isempty(idx)
        variableName = vars(idx);
    end
end

function output = normalize_name(input)
    output = lower(regexprep(string(input), '[^a-zA-Z0-9]', ''));
end

function x = as_numeric(x)

    if isnumeric(x)
        x = double(x(:));
    elseif islogical(x)
        x = double(x(:));
    elseif iscategorical(x)
        x = str2double(string(x(:)));
    elseif isstring(x) || ischar(x) || iscell(x)
        x = str2double(string(x(:)));
    else
        error('Unsupported questionnaire variable type: %s', class(x));
    end
end

function assert_variables_exist(T, variableNames)

    missing = variableNames(~ismember(variableNames, T.Properties.VariableNames));

    if ~isempty(missing)
        error('Missing required variable(s): %s', strjoin(string(missing), ', '));
    end
end

function p = exact_mcnemar_two_sided(b, c)

    validateattributes(b, {'numeric'}, {'scalar','integer','nonnegative'});
    validateattributes(c, {'numeric'}, {'scalar','integer','nonnegative'});

    nDiscordant = b + c;

    if nDiscordant == 0
        p = 1;
        return;
    end

    k = min(b, c);
    lowerTail = 0;

    for x = 0:k
        logProbability = log_nchoosek(nDiscordant, x) ...
            - nDiscordant * log(2);
        lowerTail = lowerTail + exp(logProbability);
    end

    p = min(1, 2 * lowerTail);
end

function [p, oddsRatio] = fisher_exact_two_sided(T)

    if ~isequal(size(T), [2 2])
        error('Fisher input must be a 2x2 matrix.');
    end

    if any(T(:) < 0) || any(mod(T(:),1) ~= 0)
        error('Fisher table values must be nonnegative integers.');
    end

    a = T(1,1);
    b = T(1,2);
    c = T(2,1);
    d = T(2,2);

    if b*c == 0
        if a*d > 0
            oddsRatio = Inf;
        elseif a*d == 0
            oddsRatio = NaN;
        else
            oddsRatio = 0;
        end
    else
        oddsRatio = (a*d) / (b*c);
    end

    row1 = a + b;
    col1 = a + c;
    col2 = b + d;
    total = row1 + c + d;

    xMin = max(0, row1 - col2);
    xMax = min(row1, col1);
    support = xMin:xMax;

    logDenominator = log_nchoosek(total, row1);
    probabilities = zeros(size(support));

    for i = 1:numel(support)
        x = support(i);
        probabilities(i) = exp( ...
            log_nchoosek(col1, x) + ...
            log_nchoosek(col2, row1 - x) - ...
            logDenominator);
    end

    observedIndex = find(support == a, 1, 'first');
    observedProbability = probabilities(observedIndex);

    tolerance = 1e-12;
    p = sum(probabilities(probabilities <= observedProbability + tolerance));
    p = min(max(p, 0), 1);
end

function value = log_nchoosek(n, k)

    if k < 0 || k > n
        value = -Inf;
    else
        value = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1);
    end
end

function adjusted = holm_adjust(p)

    p = p(:);
    adjusted = nan(size(p));

    valid = isfinite(p);
    pValid = p(valid);

    if isempty(pValid)
        return;
    end

    [sortedP, order] = sort(pValid);
    m = numel(sortedP);
    sortedAdjusted = zeros(m,1);

    runningMaximum = 0;

    for i = 1:m
        candidate = (m - i + 1) * sortedP(i);
        runningMaximum = max(runningMaximum, candidate);
        sortedAdjusted(i) = min(1, runningMaximum);
    end

    unsortedAdjusted = nan(m,1);
    unsortedAdjusted(order) = sortedAdjusted;
    adjusted(valid) = unsortedAdjusted;
end

function fig = make_transition_figure( ...
    N, n11, n10, n01, n00, ...
    nE1Max, nE2Max, pctE1Max, pctE2Max, ...
    differencePP, pMcNemar)

    fig = figure( ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[100 100 1150 700]);

    ax = axes(fig);
    hold(ax, 'on');

    xLeft = 1;
    xRight = 2;
    yMaximum = 1;
    yBelow = 0;

    maxCount = max([n11 n10 n01 n00 1]);

    % Calm publication-oriented palette.
    colRemainMax = [0.22 0.38 0.58];
    colDecrease  = [0.85 0.45 0.16];
    colIncrease  = [0.24 0.62 0.68];
    colRemainLow = [0.55 0.55 0.55];

    % Draw the four aggregated participant trajectories.
    draw_transition(ax, xLeft, yMaximum, xRight, yMaximum, ...
        n11, maxCount, colRemainMax);
    draw_transition(ax, xLeft, yBelow, xRight, yBelow, ...
        n00, maxCount, colRemainLow);
    draw_transition(ax, xLeft, yMaximum, xRight, yBelow, ...
        n10, maxCount, colDecrease);
    draw_transition(ax, xLeft, yBelow, xRight, yMaximum, ...
        n01, maxCount, colIncrease);

    % Category nodes.
    nodeText(ax, xLeft, yMaximum, ...
        sprintf('Maximum RAS\n%d/%d (%.1f%%)', nE1Max, N, pctE1Max));
    nodeText(ax, xLeft, yBelow, ...
        sprintf('Below maximum\n%d/%d (%.1f%%)', N-nE1Max, N, 100-pctE1Max));

    nodeText(ax, xRight, yMaximum, ...
        sprintf('Maximum RAS\n%d/%d (%.1f%%)', nE2Max, N, pctE2Max));
    nodeText(ax, xRight, yBelow, ...
        sprintf('Below maximum\n%d/%d (%.1f%%)', N-nE2Max, N, 100-pctE2Max));

    % Transition labels.
    text(ax, 1.50, 1.105, sprintf('Remained maximum: n = %d', n11), ...
        'HorizontalAlignment','center', ...
        'FontSize',12, 'FontWeight','bold', ...
        'Color',colRemainMax, 'BackgroundColor','w', 'Margin',3);

    text(ax, 1.50, -0.105, sprintf('Remained below: n = %d', n00), ...
        'HorizontalAlignment','center', ...
        'FontSize',12, 'FontWeight','bold', ...
        'Color',colRemainLow, 'BackgroundColor','w', 'Margin',3);

    text(ax, 1.42, 0.59, sprintf('Decreased: n = %d', n10), ...
        'HorizontalAlignment','center', ...
        'FontSize',12, 'FontWeight','bold', ...
        'Color',colDecrease, 'BackgroundColor','w', 'Margin',3);

    text(ax, 1.58, 0.41, sprintf('Increased: n = %d', n01), ...
        'HorizontalAlignment','center', ...
        'FontSize',12, 'FontWeight','bold', ...
        'Color',colIncrease, 'BackgroundColor','w', 'Margin',3);

    % Column titles.
    text(ax, xLeft, 1.31, 'Counting (E1)', ...
        'HorizontalAlignment','center', ...
        'FontSize',17, 'FontWeight','bold');

    text(ax, xRight, 1.31, 'Visual interference (E2)', ...
        'HorizontalAlignment','center', ...
        'FontSize',17, 'FontWeight','bold');

    % Main statistical annotation.
    resultString = sprintf([ ...
        'Maximum-awareness difference (E2 - E1): %.1f percentage points   |   ' ...
        'Two-sided exact McNemar p = %.3f'], ...
        differencePP, pMcNemar);

    text(ax, 1.50, -0.30, resultString, ...
        'HorizontalAlignment','center', ...
        'FontSize',13, 'FontWeight','bold');

    title(ax, 'Participant-level transitions in regularity awareness', ...
        'FontSize',20, 'FontWeight','bold');

    subtitle(ax, 'Maximum awareness was defined as a total RAS of 4/4', ...
        'FontSize',13);

    xlim(ax, [0.65 2.35]);
    ylim(ax, [-0.40 1.42]);
    axis(ax, 'off');
    set(ax, 'FontName', 'Arial');

    hold(ax, 'off');
end

function draw_transition(ax, x1, y1, x2, y2, count, maxCount, color)

    if count <= 0
        return;
    end

    t = linspace(0,1,250);
    smooth = t.^2 .* (3 - 2*t);

    x = x1 + (x2 - x1) .* t;
    y = y1 + (y2 - y1) .* smooth;

    minWidth = 2.0;
    maxWidth = 15.0;
    lineWidth = minWidth + (count / maxCount) * (maxWidth - minWidth);

    % White underlay separates crossing trajectories.
    plot(ax, x, y, '-', ...
        'Color','w', ...
        'LineWidth',lineWidth + 4);

    plot(ax, x, y, '-', ...
        'Color',color, ...
        'LineWidth',lineWidth);
end

function nodeText(ax, x, y, label)

    text(ax, x, y, label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',12, ...
        'FontWeight','bold', ...
        'BackgroundColor','w', ...
        'EdgeColor',[0.20 0.20 0.20], ...
        'LineWidth',1.2, ...
        'Margin',9);
end

function output = format_number(value)

    if isnan(value)
        output = 'NaN';
    elseif isinf(value)
        if value > 0
            output = 'Inf';
        else
            output = '-Inf';
        end
    else
        output = sprintf('%.3f', value);
    end
end
