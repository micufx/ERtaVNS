function OUT = eeg_export_total_scatter_4ch(outroot, grandFold, GROUPPMAP, rasFile, excludeSubs, targets, chanList)
% eeg_export_total_scatter_4ch
%
% Participant-level scatter summary for the relation between
% TOTAL regularity-awareness score and GLOBAL-EFFECT ERP mean.
%
% Main analysis:
%   - ERP global effect = Global Deviant - Global Standard
%   - Global Standard = mean(FLS, FLD)
%   - Global Deviant = mean(RLS, RLD)
%   - Parameter = mean amplitude within the GLOBAL blue-bar window
%
% Windows:
%   - E1 (Counting): uses GROUPPMAP.E1.mGlo full span
%   - E2 (Visual interference): uses GROUPPMAP.E2.mGlo full span
%   - If E2 has no global blue bar, reuses the E1 window
%
% Channels:
%   - Default representative channels: Fz, FC1, Cz, Pz
%
% Output figure:
%   - 2 x 4 participant scatter figure
%   - Top row: Counting
%   - Bottom row: Visual interference
%   - Columns: Fz, FC1, Cz, Pz
%   - Each point = one participant
%   - Participant indices are shown next to points
%
% Stats:
%   - Spearman correlation of Total vs ERP mean per channel
%   - Holm correction across the 4 channels within each experiment
%
% Tables:
%   - subject values per experiment
%   - windows table
%   - stats summary table
%
% Miguel Contreras Altamirano, 2026
% Adapted to Total-only participant-level channel scatter plots.

%% -------------------- defaults --------------------
if nargin < 5 || isempty(excludeSubs)
    excludeSubs = ["S21","S12","S23","S02","S29"];
end
if nargin < 6 || isempty(targets)
    targets = [ ...
        "Frequent Local Standard", ...
        "Rare Local Standard", ...
        "Frequent Local Deviant", ...
        "Rare Local Deviant"];
end
if nargin < 7 || isempty(chanList)
    chanList = ["Fz","FC1","Cz","Pz"];
end

excludeSubs = string(excludeSubs(:));
targets     = string(targets(:));
chanList    = string(chanList(:));

if ~exist(outroot, 'dir')
    error('outroot does not exist: %s', outroot);
end
if ~exist(grandFold, 'dir')
    mkdir(grandFold);
end
if exist('pop_loadset','file') ~= 2
    error('pop_loadset not found on MATLAB path. Make sure EEGLAB is loaded first.');
end
if nargin < 4 || isempty(rasFile) || ~exist(rasFile,'file')
    error('rasFile not found: %s', string(rasFile));
end

%% -------------------- style/config --------------------
cfg = struct();

cfg.labelE1 = "Counting";
cfg.labelE2 = "Visual interference";

cfg.colorE1 = [0 0.45 0.85];
cfg.colorE2 = [0.90 0.45 0.10];
cfg.grayLine = [0.62 0.62 0.62];

cfg.bgERP = [0.91 0.97 0.91];

cfg.dpi     = 300;
cfg.savePDF = false;

% Same width/style feel as your rainbox figures, but taller for 2x4 panels
cfg.figW = 13.8;
cfg.figH = 7.9;

cfg.pointSize = 42;
cfg.pointFaceAlpha = 0.45;
cfg.pointEdgeAlpha = 0.45;
cfg.lineWidth = 2.0;
cfg.violinWidth = 0.22;
cfg.violinAlpha = 0.16;
cfg.meanLineWidth = 3.0;
cfg.axisFontSize = 11;
cfg.titleFontSize = 12;
cfg.labelFontSize = 12;
cfg.annFontSize = 10;
cfg.idFontSize = 8;

cfg.corrType = 'Spearman';
cfg.multCompMethod = "holm";

fprintf('\n[ERP Total Scatter] Channels = %s\n', strjoin(cellstr(chanList'), ', '));
fprintf('[ERP Total Scatter] Questionnaire variable = TOTAL only\n');
fprintf('[ERP Total Scatter] Multiple-comparison correction = %s within each experiment across the %d channels\n', ...
    upper(char(cfg.multCompMethod)), numel(chanList));

%% -------------------- sanity checks --------------------
mustHaveFields = {'E1','E2'};
for k = 1:numel(mustHaveFields)
    if ~isfield(GROUPPMAP, mustHaveFields{k})
        error('GROUPPMAP.%s is missing.', mustHaveFields{k});
    end
end
req = {'times','mGlo'};
for k = 1:numel(req)
    if ~isfield(GROUPPMAP.E1, req{k})
        error('GROUPPMAP.E1.%s is missing.', req{k});
    end
    if ~isfield(GROUPPMAP.E2, req{k})
        error('GROUPPMAP.E2.%s is missing.', req{k});
    end
end

%% -------------------- windows from EEG blue bars --------------------
fallbackE1 = ex_get_fallback_window(GROUPPMAP.E1.times);
fallbackE2 = ex_get_fallback_window(GROUPPMAP.E2.times);

Wexp.Global.E1 = ex_mask_span_window(GROUPPMAP.E1.times, GROUPPMAP.E1.mGlo, fallbackE1);

% If E2 has no global blue bar, use E1 global window
if any(logical(GROUPPMAP.E2.mGlo(:)))
    Wexp.Global.E2 = ex_mask_span_window(GROUPPMAP.E2.times, GROUPPMAP.E2.mGlo, fallbackE2);
    Wexp.Global.E2.usedCountingWindow = false;
else
    Wexp.Global.E2 = Wexp.Global.E1;
    Wexp.Global.E2.usedCountingWindow = true;
end

fprintf('\n[ERP Total Scatter] WINDOWS USED FOR EXTRACTION\n');
fprintf('  Global | %s: %.0f..%.0f ms\n', ...
    cfg.labelE1, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);
fprintf('  Global | %s: %.0f..%.0f ms', ...
    cfg.labelE2, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);
if Wexp.Global.E2.usedCountingWindow
    fprintf('  (reused from Counting)\n');
else
    fprintf('\n');
end

%% -------------------- questionnaire table --------------------
Q = ex_read_ras_table(rasFile, excludeSubs);

%% -------------------- extract subject values --------------------
Te1 = ex_build_subject_table(outroot, 'E1', cfg.labelE1, Q, excludeSubs, targets, chanList, Wexp.Global.E1);
Te2 = ex_build_subject_table(outroot, 'E2', cfg.labelE2, Q, excludeSubs, targets, chanList, Wexp.Global.E2);

if isempty(Te1) && isempty(Te2)
    error('No valid subjects found after matching EEG files and questionnaire data.');
end

outSubCsvE1 = fullfile(grandFold, 'ERP_TOTAL_SCATTER_subject_values_E1.csv');
outSubCsvE2 = fullfile(grandFold, 'ERP_TOTAL_SCATTER_subject_values_E2.csv');
outSubMat   = fullfile(grandFold, 'ERP_TOTAL_SCATTER_subject_values.mat');
if ~isempty(Te1), writetable(Te1, outSubCsvE1); end
if ~isempty(Te2), writetable(Te2, outSubCsvE2); end
save(outSubMat, 'Te1', 'Te2', 'Wexp', 'cfg', '-v7.3');

%% -------------------- stats tables --------------------
TstatsE1 = ex_corr_table(Te1, chanList, cfg.labelE1, cfg.multCompMethod);
TstatsE2 = ex_corr_table(Te2, chanList, cfg.labelE2, cfg.multCompMethod);
Tstats = [TstatsE1; TstatsE2];

outStatsCsv = fullfile(grandFold, 'ERP_TOTAL_SCATTER_stats_summary.csv');
outStatsMat = fullfile(grandFold, 'ERP_TOTAL_SCATTER_stats_summary.mat');
writetable(Tstats, outStatsCsv);
save(outStatsMat, 'Tstats', 'TstatsE1', 'TstatsE2', '-v7.3');

disp(Tstats);

%% -------------------- windows table --------------------
Twin = ex_build_window_table(Wexp, cfg);
outWinCsv = fullfile(grandFold, 'ERP_TOTAL_SCATTER_sig_windows.csv');
outWinMat = fullfile(grandFold, 'ERP_TOTAL_SCATTER_sig_windows.mat');
writetable(Twin, outWinCsv);
save(outWinMat, 'Twin', 'Wexp', '-v7.3');

%% -------------------- figure --------------------
fig = figure('Color','w', 'Units','inches', 'Position',[0.5 0.5 cfg.figW cfg.figH], 'InvertHardcopy','off');
tl = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
tl.OuterPosition = [0.02 0.05 0.98 0.91];

% Common y-limits per row for comparability within experiment
row1Vals = []; row2Vals = [];
if ~isempty(Te1)
    for c = 1:numel(chanList), row1Vals = [row1Vals; Te1.(matlab.lang.makeValidName(chanList(c)))(:)]; end %#ok<AGROW>
end
if ~isempty(Te2)
    for c = 1:numel(chanList), row2Vals = [row2Vals; Te2.(matlab.lang.makeValidName(chanList(c)))(:)]; end %#ok<AGROW>
end
row1Lim = ex_nice_ylim(row1Vals);
row2Lim = ex_nice_ylim(row2Vals);

for c = 1:numel(chanList)
    ch = chanList(c);
    vname = matlab.lang.makeValidName(ch);

    % --- Counting row ---
    ax = nexttile(tl, c); hold(ax,'on');
    ex_style_axis(ax, cfg.bgERP, cfg.axisFontSize);
    if ~isempty(Te1)
        [rho,p,n] = ex_spearman(Te1.Total, Te1.(vname));
        pc = TstatsE1.CorrectedPValue(TstatsE1.Channel == ch);
        ex_plot_violin_by_total(ax, Te1.Total, Te1.(vname), cfg.colorE1, cfg);
        ex_plot_scatter(ax, Te1.Total, Te1.(vname), Te1.Index, cfg.colorE1, cfg, row1Lim);
        ex_add_fitline(ax, Te1.Total, Te1.(vname), cfg.colorE1, cfg.lineWidth);
        ex_add_statbox(ax, rho, p, pc, n, cfg.annFontSize);
    else
        title(ax, sprintf('%s', upper(char(ch))), 'FontWeight','bold');
        ylim(ax,row1Lim);
    end
    title(ax, sprintf('%s', upper(char(ch))), 'FontWeight','bold', 'FontSize', cfg.titleFontSize);
    if c == 1
        ylabel(ax, sprintf('%s\nMean global effects (\muV)', cfg.labelE1), 'FontSize', cfg.labelFontSize);
    else
        ylabel(ax, 'Global effects (\muV)', 'FontSize', cfg.labelFontSize);
    end
    xlabel(ax, 'RAS', 'FontSize', cfg.labelFontSize);

    % --- Visual interference row ---
    ax = nexttile(tl, 4 + c); hold(ax,'on');
    ex_style_axis(ax, cfg.bgERP, cfg.axisFontSize);
    if ~isempty(Te2)
        [rho,p,n] = ex_spearman(Te2.Total, Te2.(vname));
        pc = TstatsE2.CorrectedPValue(TstatsE2.Channel == ch);
        ex_plot_violin_by_total(ax, Te2.Total, Te2.(vname), cfg.colorE2, cfg);
        ex_plot_scatter(ax, Te2.Total, Te2.(vname), Te2.Index, cfg.colorE2, cfg, row2Lim);
        ex_add_fitline(ax, Te2.Total, Te2.(vname), cfg.colorE2, cfg.lineWidth);
        ex_add_statbox(ax, rho, p, pc, n, cfg.annFontSize);
    else
        title(ax, sprintf('%s', upper(char(ch))), 'FontWeight','bold');
        ylim(ax,row2Lim);
    end
    title(ax, sprintf('%s', upper(char(ch))), 'FontWeight','bold', 'FontSize', cfg.titleFontSize);
    if c == 1
        ylabel(ax, sprintf('%s\nGlobal effects (\muV)', cfg.labelE2), 'FontSize', cfg.labelFontSize);
    else
        ylabel(ax, 'Global effects (\muV)', 'FontSize', cfg.labelFontSize);
    end
    xlabel(ax, 'RAS', 'FontSize', cfg.labelFontSize);
end

sgtitle(fig, 'Global-effects (mean difference in ERP) vs Regularity Awareness Score (RAS)', ...
    'FontWeight','bold', 'FontSize', 14);

outPng = fullfile(grandFold, 'ERP_TOTAL_SCATTER_4ch_2x4.png');
exportgraphics(fig, outPng, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdf = fullfile(grandFold, 'ERP_TOTAL_SCATTER_4ch_2x4.pdf');
    exportgraphics(fig, outPdf, 'BackgroundColor', 'white');
else
    outPdf = "";
end
close(fig);

%% -------------------- outputs --------------------
OUT = struct();
OUT.subjectTableE1 = Te1;
OUT.subjectTableE2 = Te2;
OUT.windowTable    = Twin;
OUT.statsTable     = Tstats;
OUT.windows        = Wexp;
OUT.cfg            = cfg;

OUT.files.figurePNG      = outPng;
OUT.files.figurePDF      = outPdf;
OUT.files.subjectCsvE1   = outSubCsvE1;
OUT.files.subjectCsvE2   = outSubCsvE2;
OUT.files.subjectMat     = outSubMat;
OUT.files.statsCsv       = outStatsCsv;
OUT.files.statsMat       = outStatsMat;
OUT.files.windowCsv      = outWinCsv;
OUT.files.windowMat      = outWinMat;

fprintf('\n[ERP Total Scatter] Saved:\n');
fprintf('  %s\n', outPng);
fprintf('  %s\n', outSubCsvE1);
fprintf('  %s\n', outSubCsvE2);
fprintf('  %s\n', outStatsCsv);

end

%% ========================================================================
% HELPERS
% ========================================================================

function fallbackWin = ex_get_fallback_window(times)
times = double(times(:)');
if isempty(times)
    fallbackWin = [0 700];
else
    fallbackWin = [times(1) times(end)];
end
end

function W = ex_mask_span_window(relTimes, sigMask, fallbackWin)
relTimes = double(relTimes(:)');
sigMask  = logical(sigMask(:)');

W = struct();
W.hasSig = any(sigMask);
if ~W.hasSig
    W.startMS = fallbackWin(1);
    W.endMS   = fallbackWin(2);
    W.usedFallback = true;
else
    idx = find(sigMask);
    W.startMS = relTimes(idx(1));
    W.endMS   = relTimes(idx(end));
    W.usedFallback = false;
end
W.idx = relTimes >= W.startMS & relTimes <= W.endMS;
end

function T = ex_read_ras_table(rasFile, excludeSubs)
% Supports the user's current wide CSV format:
%   cols 1:6   = Counting (E1): Subject,Q1,Q2,Q3,Q4,Total
%   cols 10:15 = Visual interference (E2): Subject,Q1,Q2,Q3,Q4,Total
% Falls back to a tidy format if those columns are not present.

raw = readcell(rasFile, 'FileType','text');
if size(raw,1) >= 3 && size(raw,2) >= 15
    hdr1 = string(raw(1,:));
    hdr2 = string(raw(2,:));

    hasWideE1 = any(contains(lower(hdr1), "active counting")) || (strtrim(lower(hdr2(1))) == "subject");
    hasWideE2 = any(contains(lower(hdr1), "dual task")) || (strtrim(lower(hdr2(10))) == "subject");

    if hasWideE1 && hasWideE2
        data = raw(3:end,:);

        subjE1  = ex_normalize_sub_ids(string(data(:,1)));
        q1E1    = ex_to_numeric(data(:,2));
        q2E1    = ex_to_numeric(data(:,3));
        q3E1    = ex_to_numeric(data(:,4));
        q4E1    = ex_to_numeric(data(:,5));
        totE1   = ex_to_numeric(data(:,6));

        subjE2  = ex_normalize_sub_ids(string(data(:,10)));
        q1E2    = ex_to_numeric(data(:,11));
        q2E2    = ex_to_numeric(data(:,12));
        q3E2    = ex_to_numeric(data(:,13));
        q4E2    = ex_to_numeric(data(:,14));
        totE2   = ex_to_numeric(data(:,15));

        T1 = table(subjE1, repmat("E1", numel(subjE1), 1), q1E1, q2E1, q3E1, q4E1, totE1, ...
            'VariableNames', {'Subject','Experiment','Q1','Q2','Q3','Q4','Total'});
        T2 = table(subjE2, repmat("E2", numel(subjE2), 1), q1E2, q2E2, q3E2, q4E2, totE2, ...
            'VariableNames', {'Subject','Experiment','Q1','Q2','Q3','Q4','Total'});

        T = [T1; T2];
        T = T(~ismember(lower(T.Subject), lower(excludeSubs)), :);
        T = T(~ismissing(T.Subject) & strlength(strtrim(T.Subject)) > 0, :);
        validTotal = isfinite(T.Total) | isfinite(T.Q1) | isfinite(T.Q2) | isfinite(T.Q3) | isfinite(T.Q4);
        T = T(validTotal, :);

        missTot = ~isfinite(T.Total);
        if any(missTot)
            T.Total(missTot) = T.Q1(missTot) + T.Q2(missTot) + T.Q3(missTot) + T.Q4(missTot);
        end

        [~,ia] = unique(T(:,{'Subject','Experiment'}), 'rows', 'stable');
        T = T(ia,:);
        return;
    end
end

% ---------- fallback tidy parser ----------
T0 = readtable(rasFile, 'VariableNamingRule','preserve');
vars = string(T0.Properties.VariableNames);

subCol = ex_find_first_match(vars, ["subject","participant","sub","id"]);
q1Col  = ex_find_qcol(vars, 1);
q2Col  = ex_find_qcol(vars, 2);
q3Col  = ex_find_qcol(vars, 3);
q4Col  = ex_find_qcol(vars, 4);
totCol = ex_find_total_col(vars);
expCol = ex_find_first_match(vars, ["experiment","exp","task","group","condition"]);

if isempty(subCol)
    error('Could not find subject column in rasFile. Variables are: %s', strjoin(cellstr(vars), ', '));
end
if isempty(q1Col) || isempty(q2Col) || isempty(q3Col) || isempty(q4Col)
    error('Could not find Q1-Q4 columns in rasFile. Variables are: %s', strjoin(cellstr(vars), ', '));
end

Subject = ex_normalize_sub_ids(string(T0.(subCol)));
Q1 = ex_to_numeric(T0.(q1Col));
Q2 = ex_to_numeric(T0.(q2Col));
Q3 = ex_to_numeric(T0.(q3Col));
Q4 = ex_to_numeric(T0.(q4Col));

if isempty(totCol)
    Total = Q1 + Q2 + Q3 + Q4;
else
    Total = ex_to_numeric(T0.(totCol));
end

if isempty(expCol)
    Experiment = repmat("E1", numel(Subject), 1);
else
    Experiment = string(T0.(expCol));
    Experiment = regexprep(lower(strtrim(Experiment)), '\s+', '');
    Experiment(startsWith(Experiment, "e1") | contains(Experiment,"count")) = "E1";
    Experiment(startsWith(Experiment, "e2") | contains(Experiment,"dual") | contains(Experiment,"visual")) = "E2";
end

T = table(Subject, Experiment, Q1, Q2, Q3, Q4, Total);
T = T(~ismember(lower(T.Subject), lower(excludeSubs)), :);
T = T(~ismissing(T.Subject) & strlength(T.Subject) > 0, :);
[~,ia] = unique(T(:,{'Subject','Experiment'}), 'rows', 'stable');
T = T(ia,:);
end

function name = ex_find_first_match(vars, patterns)
name = string.empty(0,1);
varsLow = lower(vars);
for p = 1:numel(patterns)
    idx = find(contains(varsLow, lower(patterns(p))), 1, 'first');
    if ~isempty(idx)
        name = string(vars(idx));
        return;
    end
end
end

function name = ex_find_qcol(vars, qnum)
name = string.empty(0,1);
varsLow = lower(vars);
patterns = [ ...
    "q" + qnum, ...
    "q_" + qnum, ...
    "q " + qnum, ...
    "question" + qnum, ...
    "question_" + qnum, ...
    "question " + qnum];
for i = 1:numel(vars)
    v = varsLow(i);
    ok = false;
    for p = 1:numel(patterns)
        if contains(v, patterns(p))
            ok = true; break;
        end
    end
    if ok
        name = vars(i);
        return;
    end
end
end

function name = ex_find_total_col(vars)
name = string.empty(0,1);
varsLow = lower(vars);
patterns = ["total","sum","regularity awareness score","ras"];
for i = 1:numel(vars)
    v = varsLow(i);
    for p = 1:numel(patterns)
        if contains(v, patterns(p)) && ~contains(v,'q')
            name = vars(i);
            return;
        end
    end
end
end

function S = ex_normalize_sub_ids(x)
x = string(x(:));
S = strings(size(x));
for i = 1:numel(x)
    xi = strtrim(x(i));
    tok = regexp(char(xi), '(\d+)', 'tokens', 'once');
    if isempty(tok)
        S(i) = "";
    else
        S(i) = "S" + pad(string(str2double(tok{1})), 2, 'left', '0');
    end
end
end

function y = ex_to_numeric(x)
if isnumeric(x)
    y = double(x(:));
elseif islogical(x)
    y = double(x(:));
elseif iscell(x)
    y = nan(numel(x),1);
    for i = 1:numel(x)
        y(i) = str2double(string(x{i}));
    end
else
    y = str2double(string(x(:)));
end
end

function Texp = ex_build_subject_table(outroot, expKey, expLabel, Q, excludeSubs, targets, chanList, W)
subs = ex_list_subjects_with_file(outroot, expKey, excludeSubs);
rows = cell(0,1);

for i = 1:numel(subs)
    sid = subs(i);
    iq = find(Q.Subject == sid & Q.Experiment == string(expKey), 1, 'first');
    if isempty(iq), continue; end

    setFile = fullfile(outroot, char(sid), sprintf('auditory_all_LG_%s_%s.set', char(sid), expKey));
    if ~exist(setFile,'file'), continue; end

    vals = ex_load_subject_globaleffect_means(setFile, chanList, targets, W.startMS, W.endMS);
    if isempty(vals), continue; end

    row = struct();
    row.Subject = sid;
    row.Q1 = Q.Q1(iq);
    row.Q2 = Q.Q2(iq);
    row.Q3 = Q.Q3(iq);
    row.Q4 = Q.Q4(iq);
    row.Total = Q.Total(iq);
    for c = 1:numel(chanList)
        vname = matlab.lang.makeValidName(chanList(c));
        row.(vname) = vals.(vname);
    end
    rows{end+1,1} = row; %#ok<AGROW>
end

if isempty(rows)
    Texp = table();
    return;
end

Texp = struct2table([rows{:}]');
Texp = sortrows(Texp, 'Subject');
Texp.Index = (1:height(Texp))';
Texp = movevars(Texp, 'Index', 'Before', 'Subject');
Texp.Experiment = repmat(string(expLabel), height(Texp), 1);
Texp = movevars(Texp, 'Experiment', 'After', 'Subject');
end

function subs = ex_list_subjects_with_file(outroot, expKey, excludeSubs)
d = dir(fullfile(outroot, 'S*'));
d = d([d.isdir]);
subs = string({d.name});
subs = subs(~ismember(subs, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(subs), '^S\d+$', 'once'));
subs = subs(keep);
subs = subs(~ismember(lower(subs), lower(excludeSubs)));

hasFile = false(size(subs));
for i = 1:numel(subs)
    f = fullfile(outroot, char(subs(i)), sprintf('auditory_all_LG_%s_%s.set', char(subs(i)), expKey));
    hasFile(i) = exist(f,'file') == 2;
end
subs = subs(hasFile);
subs = ex_sort_subject_ids(subs);
end

function subs = ex_sort_subject_ids(subs)
subs = string(subs);
nums = nan(size(subs));
for i = 1:numel(subs)
    tok = regexp(char(subs(i)), '^S(\d+)$', 'tokens', 'once');
    if ~isempty(tok), nums(i) = str2double(tok{1}); end
end
[~,ix] = sort(nums);
subs = subs(ix);
end

function vals = ex_load_subject_globaleffect_means(setFile, chanList, targets, startMS, endMS)
vals = [];
EEG = pop_loadset(setFile);
EEG = eeg_checkset(EEG);

if ~(isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type'))
    return;
end

M = false(4, EEG.trials);
for ep = 1:EEG.trials
    lablist = ex_to_strlist(EEG.epoch(ep).eventlg_type);
    for k = 1:4
        M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
    end
end

idxFLS = find(M(1,:));
idxRLS = find(M(2,:));
idxFLD = find(M(3,:));
idxRLD = find(M(4,:));
if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
    return;
end

vals = struct();
t = double(EEG.times(:)');
winIdx = t >= startMS & t <= endMS;
if ~any(winIdx)
    return;
end

for c = 1:numel(chanList)
    ch = find(strcmpi({EEG.chanlocs.labels}, char(chanList(c))), 1);
    if isempty(ch)
        vals.(matlab.lang.makeValidName(chanList(c))) = NaN;
        continue;
    end
    FLS = squeeze(mean(double(EEG.data(ch,:,idxFLS)),3)); FLS = FLS(:)';
    RLS = squeeze(mean(double(EEG.data(ch,:,idxRLS)),3)); RLS = RLS(:)';
    FLD = squeeze(mean(double(EEG.data(ch,:,idxFLD)),3)); FLD = FLD(:)';
    RLD = squeeze(mean(double(EEG.data(ch,:,idxRLD)),3)); RLD = RLD(:)';

    GS = mean([FLS; FLD], 1, 'omitnan');
    GD = mean([RLS; RLD], 1, 'omitnan');
    GE = GD - GS;

    vals.(matlab.lang.makeValidName(chanList(c))) = mean(GE(winIdx), 'omitnan');
end
end

function S = ex_to_strlist(v)
if isstring(v)
    S = v;
elseif ischar(v)
    S = string(v);
elseif iscell(v)
    S = strings(0,1);
    for i = 1:numel(v)
        if isstring(v{i}) || ischar(v{i})
            S(end+1,1) = string(v{i}); %#ok<AGROW>
        end
    end
else
    S = strings(0,1);
end
S = strtrim(regexprep(S,'\s+',' '));
end

function T = ex_corr_table(Texp, chanList, expLabel, multCompMethod)
if isempty(Texp)
    T = table();
    return;
end

Experiment = repmat(string(expLabel), numel(chanList), 1);
Channel = chanList(:);
N = nan(numel(chanList),1);
Rho = nan(numel(chanList),1);
PValue = nan(numel(chanList),1);
CorrectedPValue = nan(numel(chanList),1);
Significance = strings(numel(chanList),1);

for c = 1:numel(chanList)
    vname = matlab.lang.makeValidName(chanList(c));
    [rho,p,n] = ex_spearman(Texp.Total, Texp.(vname));
    Rho(c) = rho; PValue(c) = p; N(c) = n;
end
CorrectedPValue = ex_correct_pvalues(PValue, multCompMethod);
for c = 1:numel(chanList)
    Significance(c) = ex_p_to_stars(CorrectedPValue(c));
end

T = table(Experiment, Channel, N, Rho, PValue, CorrectedPValue, Significance);
end

function [rho,p,n] = ex_spearman(x,y)
m = isfinite(x) & isfinite(y);
x = x(m); y = y(m); n = numel(x);
if n < 3 || numel(unique(x)) < 2 || numel(unique(y)) < 2
    rho = NaN; p = NaN; return;
end
[R,P] = corr(x, y, 'Type','Spearman', 'Rows','complete');
rho = R; p = P;
end

function Twin = ex_build_window_table(Wexp, cfg)
Effect = ["Global"; "Global"];
Experiment = [cfg.labelE1; cfg.labelE2];
StartMS = [Wexp.Global.E1.startMS; Wexp.Global.E2.startMS];
EndMS   = [Wexp.Global.E1.endMS; Wexp.Global.E2.endMS];
UsedFallback = [Wexp.Global.E1.usedFallback; Wexp.Global.E2.usedFallback];
UsedCountingWindow = [false; Wexp.Global.E2.usedCountingWindow];
WindowType = repmat("Experiment-specific global blue bar", 2, 1);
Twin = table(Effect, WindowType, Experiment, StartMS, EndMS, UsedFallback, UsedCountingWindow);
end

function ex_style_axis(ax, bgColor, fontSize)
set(ax, 'Color', bgColor, 'FontSize', fontSize, 'LineWidth', 1, ...
    'Box','off', 'Layer','top', 'TickDir','out');
grid(ax, 'on');
end

function ex_plot_scatter(ax, x, y, ids, rgb, cfg, yLim)
m = isfinite(x) & isfinite(y);
x = x(m); y = y(m); ids = ids(m);
if isempty(x)
    ylim(ax, yLim);
    return;
end

rng(1);
xj = x + (rand(size(x)) - 0.5) * 0.06;
scatter(ax, xj, y, cfg.pointSize, ...
    'MarkerFaceColor', rgb, ...
    'MarkerEdgeColor', [0.30 0.30 0.30], ...
    'LineWidth', 0.35, ...
    'MarkerFaceAlpha', cfg.pointFaceAlpha, ...
    'MarkerEdgeAlpha', cfg.pointEdgeAlpha, ...
    'HandleVisibility', 'off');

xRange = range(xlim(ax));
if xRange == 0, xRange = 1; end
% for i = 1:numel(x)
%     text(ax, xj(i) + 0.02*xRange, y(i), sprintf('%d', ids(i)), ...
%         'FontSize', cfg.idFontSize, 'Color', [0.20 0.20 0.20], ...
%         'VerticalAlignment','middle', 'Clipping','on');
% end

xt = unique(x(~isnan(x)));
if isempty(xt), xt = unique(round(x)); end
set(ax, 'XTick', xt);
ylim(ax, yLim);
end

function ex_plot_violin_by_total(ax, x, y, rgb, cfg)
% Full violin per Total category, using the same experiment color.
m = isfinite(x) & isfinite(y);
x = x(m); y = y(m);
if isempty(x), return; end

ux = unique(x(:)');
ux = ux(isfinite(ux));
for k = 1:numel(ux)
    xc = ux(k);
    yy = y(x == xc);
    yy = yy(isfinite(yy));
    if numel(yy) < 2
        if numel(yy) == 1
            plot(ax, [xc-cfg.violinWidth/2, xc+cfg.violinWidth/2], [yy yy], '-', ...
                'Color', rgb, 'LineWidth', cfg.meanLineWidth, 'HandleVisibility','off');
        end
        continue;
    end
    try
        if numel(unique(yy)) == 1
            ygrid = linspace(yy(1)-eps(yy(1)+1), yy(1)+eps(yy(1)+1), 64);
            f = ones(size(ygrid));
        else
            [f, ygrid] = ksdensity(yy, 'NumPoints', 200);
        end
    catch
        [N, edges] = histcounts(yy, 30, 'Normalization', 'pdf');
        ygrid = edges(1:end-1) + diff(edges)/2;
        f = N;
    end
    if isempty(f) || all(~isfinite(f)) || max(f) <= 0
        continue;
    end
    f = f ./ max(f) * cfg.violinWidth;
    xpoly = [xc - f, fliplr(xc + f)];
    ypoly = [ygrid, fliplr(ygrid)];
    patch(ax, xpoly, ypoly, rgb, ...
        'FaceAlpha', cfg.violinAlpha, ...
        'EdgeColor', rgb, ...
        'LineWidth', 1.4, ...
        'HandleVisibility','off');
    mu = mean(yy, 'omitnan');
    plot(ax, [xc-cfg.violinWidth/2, xc+cfg.violinWidth/2], [mu mu], '-', ...
        'Color', rgb, 'LineWidth', cfg.meanLineWidth, 'HandleVisibility','off');
end
end

function ex_add_fitline(ax, x, y, rgb, lw)
m = isfinite(x) & isfinite(y);
x = x(m); y = y(m);
if numel(x) < 2 || numel(unique(x)) < 2
    return;
end
p = polyfit(x, y, 1);
xf = linspace(min(x), max(x), 100);
yf = polyval(p, xf);
plot(ax, xf, yf, '-', 'Color', rgb, 'LineWidth', lw, 'HandleVisibility','off');
end

function ex_add_statbox(ax, rho, p, pc, n, fontSize)

if ~isfinite(rho)
    s1 = '\rho = n/a';
else
    s1 = sprintf('\\rho = %.3f', rho);
end

if ~isfinite(p)
    s2 = 'p = n/a';
else
    s2 = sprintf('p = %.3f', p);
end

if ~isfinite(pc)
    s3 = 'p_c = n/a';
else
    s3 = sprintf('p_c = %.3f', pc);
end

s4 = sprintf('n = %d', n);

txt = sprintf('%s\n%s\n%s\n%s', s1, s2, s3, s4);

text(ax, 0.03, 0.97, txt, ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize', fontSize, ...
    'FontWeight','bold', ...
    'BackgroundColor',[1 1 1], ...
    'Margin',4, ...
    'EdgeColor',[0.85 0.85 0.85]);

end
function yLim = ex_nice_ylim(yData)
yData = yData(isfinite(yData));
if isempty(yData)
    yLim = [-1 1];
    return;
end
yMin = min(yData); yMax = max(yData); yr = yMax - yMin;
if yr <= 0, yr = max(abs(yMax), 1); end
padLo = 0.12 * yr; padHi = 0.16 * yr;
yLim = [yMin - padLo, yMax + padHi];
end

function txt = ex_p_to_stars(p)
if ~isfinite(p)
    txt = "n/a";
elseif p < 0.001
    txt = "***";
elseif p < 0.01
    txt = "**";
elseif p < 0.05
    txt = "*";
else
    txt = "ns";
end
end

function pAdj = ex_correct_pvalues(pRaw, methodStr)
pAdj = nan(size(pRaw));
ok = isfinite(pRaw);
if ~any(ok), return; end
p = pRaw(ok);
switch lower(char(methodStr))
    case 'none'
        pAdj(ok) = p;
    case 'bonferroni'
        pAdj(ok) = min(p * numel(p), 1);
    otherwise
        pAdj(ok) = ex_holm_adjust(p);
end
end

function pAdj = ex_holm_adjust(p)
origSize = size(p);
p = p(:);
m = numel(p);
[ps, ord] = sort(p, 'ascend');
adjSorted = (m - (1:m)' + 1) .* ps;
adjSorted = min(adjSorted, 1);
adjSorted = cummax(adjSorted);
pAdj = nan(size(p));
pAdj(ord) = adjSorted;
pAdj = reshape(pAdj, origSize);
end
