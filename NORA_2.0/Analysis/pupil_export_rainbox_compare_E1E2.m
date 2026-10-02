function OUT = pupil_export_rainbox_compare_E1E2(outpath, grandOutFolder, STATS, excludeSubs, useField)
% pupil_export_rainbox_compare_E1E2
%
% Creates 4 publication-style figures:
%   1) PUPIL_RAINBOX_COMPARE_Local_E1vsE2.png
%   2) PUPIL_RAINBOX_COMPARE_Global_E1vsE2.png
%   3) PUPIL_RAINBOX_COMPARE_LocalSlope_E1vsE2.png
%   4) PUPIL_RAINBOX_COMPARE_GlobalSlope_E1vsE2.png
%
% Figures 1-2:
%   LEFT  = peak latency
%   RIGHT = mean pupil size (z-score)
%
% Figures 3-4:
%   LEFT  = peak latency
%   RIGHT = slope (z-score/s)
%
% Within each panel, Counting vs Visual interference are compared for:
%   Local figure  -> Local Standard, Local Deviant
%   Global figure -> Global Standard, Global Deviant
%
% Also adds the long top bracket:
%   comparison of the within-experiment Deviant-Standard contrast
%   between Counting and Visual interference.
%
% Window logic:
%   - Values are extracted using the FULL blue-bar span of EACH experiment/effect.
%   - Same experiment-specific effect window is used for both conditions of that effect.
%
% Slope logic:
%   - Default = linear slope from window START to the subject-specific PEAK
%     within that same blue-bar window.
%   - Units = z-score/s
%
% Outputs:
%   - subject values table
%   - blue-bar windows table
%   - stats summary table
%
% Miguel Contreras Altamirano, 2026

%% -------------------- defaults --------------------
if nargin < 4 || isempty(excludeSubs)
    excludeSubs = strings(0,1);
end
if nargin < 5 || isempty(useField)
    useField = "PUPIL_Z";
end

useField    = string(useField);
excludeSubs = string(excludeSubs);

if ~exist(outpath, 'dir')
    error('outpath does not exist: %s', outpath);
end
if ~exist(grandOutFolder, 'dir')
    mkdir(grandOutFolder);
end

%% -------------------- style/config --------------------
cfg = struct();

cfg.labelE1 = "Counting";
cfg.labelE2 = "Visual interference";

% SAME COLORS AS YOUR OTHER PLOTS
cfg.colorE1 = [0 0.45 0.85];
cfg.colorE2 = [0.90 0.45 0.10];

cfg.grayLine = [0.62 0.62 0.62];

% cfg.bgLatency = [0.93 0.93 0.93]; % background color left panel
% cfg.bgPupil   = [0.99 0.94 0.89]; % background color right panel

cfg.bgLatency = [1 1 1];
cfg.bgPupil   = [1 1 1];

cfg.statsMethod = "signrank";   % "signrank" | "ttest"
cfg.dpi         = 300;
cfg.savePDF     = false;

cfg.figW = 13.8;  % inches
cfg.figH = 4.8;   % inches

cfg.rainWidth     = 0.13;
cfg.boxWidth      = 0.18;
cfg.jitter        = 0.022;
cfg.pointSize     = 28;
cfg.boxLineWidth  = 2.0;
cfg.meanLineWidth = 3.0;
cfg.violinAlpha   = 0.22;

cfg.legendLocation = 'southwest';

% -------- multiple-comparison correction --------
cfg.multCompMethod = "holm";   % "holm" | "bonferroni" | "none"

fprintf(['\n[RainBox] MULTIPLE-COMPARISON CORRECTION = %s within each ' ...
    'Effect x Measure family (3 tests: Standard, Deviant, Deviant-Standard contrast)\n'], ...
    upper(char(cfg.multCompMethod)));

% -------- slope settings --------
% "start_to_peak" is recommended here because it captures dilation build-up
% without later constriction flattening the estimate.
cfg.slopeMode  = "start_to_peak";   % "start_to_peak" | "full_window"
cfg.slopeUnit  = "ms";               % "s" or "ms"

if strcmpi(cfg.slopeUnit, "ms")
    cfg.slopeLabel = "Pupil slope (z-score/ms)";
else
    cfg.slopeLabel = "Pupil slope (z-score/s)";
end

fprintf('\n[RainBox] SLOPE DEFINITION = linear slope from window start to subject peak inside that window (%s)\n', cfg.slopeUnit);

%% -------------------- sanity checks --------------------
mustHaveFields = {'E1','E2'};
for k = 1:numel(mustHaveFields)
    if ~isfield(STATS, mustHaveFields{k})
        error('STATS.%s is missing.', mustHaveFields{k});
    end
end

req = {'relTimes','maskLocal','maskGlobal'};
for k = 1:numel(req)
    if ~isfield(STATS.E1, req{k})
        error('STATS.E1.%s is missing.', req{k});
    end
    if ~isfield(STATS.E2, req{k})
        error('STATS.E2.%s is missing.', req{k});
    end
end

%% -------------------- windows from blue bars --------------------
fallbackE1 = rbx_get_fallback_window(STATS.E1);
fallbackE2 = rbx_get_fallback_window(STATS.E2);

% FULL blue-bar span per experiment/effect
Wexp.Local.E1  = rbx_mask_span_window(STATS.E1.relTimes, STATS.E1.maskLocal,  fallbackE1);
Wexp.Local.E2  = rbx_mask_span_window(STATS.E2.relTimes, STATS.E2.maskLocal,  fallbackE2);

Wexp.Global.E1 = rbx_mask_span_window(STATS.E1.relTimes, STATS.E1.maskGlobal, fallbackE1);
Wexp.Global.E2 = rbx_mask_span_window(STATS.E2.relTimes, STATS.E2.maskGlobal, fallbackE2);

fprintf('\n[RainBox] WINDOWS USED FOR EXTRACTION\n');
fprintf('  Local  | %s: %.0f..%.0f ms | %s: %.0f..%.0f ms\n', ...
    cfg.labelE1, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, ...
    cfg.labelE2, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);

fprintf('  Global | %s: %.0f..%.0f ms | %s: %.0f..%.0f ms\n', ...
    cfg.labelE1, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, ...
    cfg.labelE2, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);

%% -------------------- subjects present in BOTH E1 and E2 --------------------
subsE1 = rbx_list_subjects_with_file(outpath, 'E1', excludeSubs);
subsE2 = rbx_list_subjects_with_file(outpath, 'E2', excludeSubs);

subNames = subsE1(ismember(subsE1, subsE2));
subNames = rbx_sort_subject_ids(subNames);

if isempty(subNames)
    error('No common subjects found with both E1 and E2 pupil files.');
end

%% -------------------- subject-level extraction --------------------
nSub = numel(subNames);

% Local latency
LS_Lat_E1 = nan(nSub,1); LS_Lat_E2 = nan(nSub,1);
LD_Lat_E1 = nan(nSub,1); LD_Lat_E2 = nan(nSub,1);

% Local mean
LS_Mean_E1 = nan(nSub,1); LS_Mean_E2 = nan(nSub,1);
LD_Mean_E1 = nan(nSub,1); LD_Mean_E2 = nan(nSub,1);

% Local slope
LS_Slope_E1 = nan(nSub,1); LS_Slope_E2 = nan(nSub,1);
LD_Slope_E1 = nan(nSub,1); LD_Slope_E2 = nan(nSub,1);

% Global latency
GS_Lat_E1 = nan(nSub,1); GS_Lat_E2 = nan(nSub,1);
GD_Lat_E1 = nan(nSub,1); GD_Lat_E2 = nan(nSub,1);

% Global mean
GS_Mean_E1 = nan(nSub,1); GS_Mean_E2 = nan(nSub,1);
GD_Mean_E1 = nan(nSub,1); GD_Mean_E2 = nan(nSub,1);

% Global slope
GS_Slope_E1 = nan(nSub,1); GS_Slope_E2 = nan(nSub,1);
GD_Slope_E1 = nan(nSub,1); GD_Slope_E2 = nan(nSub,1);

for i = 1:nSub
    sID = subNames(i);

    fileE1 = fullfile(outpath, char(sID), sprintf('%s_E1_pupil_preproc.mat', sID));
    fileE2 = fullfile(outpath, char(sID), sprintf('%s_E2_pupil_preproc.mat', sID));

    if ~exist(fileE1,'file') || ~exist(fileE2,'file')
        continue;
    end

    C1 = rbx_load_subject_curves(fileE1, useField);
    C2 = rbx_load_subject_curves(fileE2, useField);

    if isempty(C1) || isempty(C2)
        continue;
    end

    % ---------- LOCAL ----------
    LS_Mean_E1(i) = rbx_mean_in_window(C1.relTimes, C1.LS, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS);
    LD_Mean_E1(i) = rbx_mean_in_window(C1.relTimes, C1.LD, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS);

    LS_Mean_E2(i) = rbx_mean_in_window(C2.relTimes, C2.LS, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);
    LD_Mean_E2(i) = rbx_mean_in_window(C2.relTimes, C2.LD, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);

    LS_Lat_E1(i) = rbx_curve_peak_latency(C1.relTimes, C1.LS, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS);
    LD_Lat_E1(i) = rbx_curve_peak_latency(C1.relTimes, C1.LD, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS);

    LS_Lat_E2(i) = rbx_curve_peak_latency(C2.relTimes, C2.LS, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);
    LD_Lat_E2(i) = rbx_curve_peak_latency(C2.relTimes, C2.LD, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);

    LS_Slope_E1(i) = rbx_curve_slope(C1.relTimes, C1.LS, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, cfg.slopeMode, cfg.slopeUnit);
    LD_Slope_E1(i) = rbx_curve_slope(C1.relTimes, C1.LD, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, cfg.slopeMode, cfg.slopeUnit);

    LS_Slope_E2(i) = rbx_curve_slope(C2.relTimes, C2.LS, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS, cfg.slopeMode, cfg.slopeUnit);
    LD_Slope_E2(i) = rbx_curve_slope(C2.relTimes, C2.LD, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS, cfg.slopeMode, cfg.slopeUnit);

    % ---------- GLOBAL ----------
    GS_Mean_E1(i) = rbx_mean_in_window(C1.relTimes, C1.GS, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);
    GD_Mean_E1(i) = rbx_mean_in_window(C1.relTimes, C1.GD, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);

    GS_Mean_E2(i) = rbx_mean_in_window(C2.relTimes, C2.GS, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);
    GD_Mean_E2(i) = rbx_mean_in_window(C2.relTimes, C2.GD, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);

    GS_Lat_E1(i) = rbx_curve_peak_latency(C1.relTimes, C1.GS, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);
    GD_Lat_E1(i) = rbx_curve_peak_latency(C1.relTimes, C1.GD, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);

    GS_Lat_E2(i) = rbx_curve_peak_latency(C2.relTimes, C2.GS, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);
    GD_Lat_E2(i) = rbx_curve_peak_latency(C2.relTimes, C2.GD, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);

    GS_Slope_E1(i) = rbx_curve_slope(C1.relTimes, C1.GS, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, cfg.slopeMode, cfg.slopeUnit);
    GD_Slope_E1(i) = rbx_curve_slope(C1.relTimes, C1.GD, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, cfg.slopeMode, cfg.slopeUnit);

    GS_Slope_E2(i) = rbx_curve_slope(C2.relTimes, C2.GS, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS, cfg.slopeMode, cfg.slopeUnit);
    GD_Slope_E2(i) = rbx_curve_slope(C2.relTimes, C2.GD, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS, cfg.slopeMode, cfg.slopeUnit);
end

%% -------------------- subject value table --------------------
Tsub = table( ...
    subNames(:), ...
    LS_Lat_E1, LS_Lat_E2, ...
    LD_Lat_E1, LD_Lat_E2, ...
    LS_Mean_E1, LS_Mean_E2, ...
    LD_Mean_E1, LD_Mean_E2, ...
    LS_Slope_E1, LS_Slope_E2, ...
    LD_Slope_E1, LD_Slope_E2, ...
    GS_Lat_E1, GS_Lat_E2, ...
    GD_Lat_E1, GD_Lat_E2, ...
    GS_Mean_E1, GS_Mean_E2, ...
    GD_Mean_E1, GD_Mean_E2, ...
    GS_Slope_E1, GS_Slope_E2, ...
    GD_Slope_E1, GD_Slope_E2, ...
    'VariableNames', { ...
    'Subject', ...
    'LocalStandard_Latency_E1_ms', 'LocalStandard_Latency_E2_ms', ...
    'LocalDeviant_Latency_E1_ms',  'LocalDeviant_Latency_E2_ms', ...
    'LocalStandard_Mean_E1',       'LocalStandard_Mean_E2', ...
    'LocalDeviant_Mean_E1',        'LocalDeviant_Mean_E2', ...
    'LocalStandard_Slope_E1',      'LocalStandard_Slope_E2', ...
    'LocalDeviant_Slope_E1',       'LocalDeviant_Slope_E2', ...
    'GlobalStandard_Latency_E1_ms','GlobalStandard_Latency_E2_ms', ...
    'GlobalDeviant_Latency_E1_ms', 'GlobalDeviant_Latency_E2_ms', ...
    'GlobalStandard_Mean_E1',      'GlobalStandard_Mean_E2', ...
    'GlobalDeviant_Mean_E1',       'GlobalDeviant_Mean_E2', ...
    'GlobalStandard_Slope_E1',     'GlobalStandard_Slope_E2', ...
    'GlobalDeviant_Slope_E1',      'GlobalDeviant_Slope_E2'});

outCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_subject_values.csv');
outMat = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_subject_values.mat');
writetable(Tsub, outCsv);
save(outMat, 'Tsub', 'Wexp', 'cfg', '-v7.3');

%% -------------------- blue-bar windows table --------------------
Twin = rbx_build_window_table(Wexp, cfg);

outWinCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_sig_windows.csv');
outWinMat = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_sig_windows.mat');
writetable(Twin, outWinCsv);
save(outWinMat, 'Twin', 'Wexp', '-v7.3');

disp(Twin);

%% -------------------- stats summary table --------------------
Tstats = rbx_build_stats_table(cfg, ...
    LS_Lat_E1, LS_Lat_E2, LD_Lat_E1, LD_Lat_E2, ...
    LS_Mean_E1, LS_Mean_E2, LD_Mean_E1, LD_Mean_E2, ...
    LS_Slope_E1, LS_Slope_E2, LD_Slope_E1, LD_Slope_E2, ...
    GS_Lat_E1, GS_Lat_E2, GD_Lat_E1, GD_Lat_E2, ...
    GS_Mean_E1, GS_Mean_E2, GD_Mean_E1, GD_Mean_E2, ...
    GS_Slope_E1, GS_Slope_E2, GD_Slope_E1, GD_Slope_E2);

outStatsCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_stats_summary.csv');
outStatsMat = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_stats_summary.mat');
writetable(Tstats, outStatsCsv);
save(outStatsMat, 'Tstats', '-v7.3');

disp(Tstats);

% ---------- corrected bracket labels for figures ----------
B_Local_Lat   = rbx_get_bracket_info(Tstats, "Local",  "Peak latency");
B_Local_Mean  = rbx_get_bracket_info(Tstats, "Local",  "Mean pupil size");
B_Local_Slope = rbx_get_bracket_info(Tstats, "Local",  "Slope pupil size");

B_Global_Lat   = rbx_get_bracket_info(Tstats, "Global", "Peak latency");
B_Global_Mean  = rbx_get_bracket_info(Tstats, "Global", "Mean pupil size");
B_Global_Slope = rbx_get_bracket_info(Tstats, "Global", "Slope pupil size");

fprintf(['\n[RainBox] COUNTING vs VISUAL INTERFERENCE SUMMARY ' ...
    '(%s-corrected within each Effect x Measure 3-test family)\n'], ...
    upper(char(cfg.multCompMethod)));

for r = 1:height(Tstats)
    fprintf('  %s | %s | %s: n=%d, raw p=%.4f, corrected p=%.4f, %s\n', ...
        Tstats.Effect(r), Tstats.Measure(r), Tstats.Condition(r), ...
        Tstats.Npairs(r), Tstats.RawPValue(r), Tstats.CorrectedPValue(r), ...
        Tstats.Significance(r));
end

%% -------------------- figures: MEAN --------------------
figLocal = rbx_make_effect_figure( ...
    "Local", ...
    ["Local Standard","Local Deviant"], ...
    {LS_Lat_E1, LD_Lat_E1}, {LS_Lat_E2, LD_Lat_E2}, ...
    {LS_Mean_E1, LD_Mean_E1}, {LS_Mean_E2, LD_Mean_E2}, ...
    "mean", "Pupil mean (z-score)", ...
    B_Local_Lat, B_Local_Mean, cfg);

outPngLocal = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_Local_E1vsE2.png');
exportgraphics(figLocal, outPngLocal, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfLocal = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_Local_E1vsE2.pdf');
    exportgraphics(figLocal, outPdfLocal, 'BackgroundColor', 'white');
end
close(figLocal);

figGlobal = rbx_make_effect_figure( ...
    "Global", ...
    ["Global Standard","Global Deviant"], ...
    {GS_Lat_E1, GD_Lat_E1}, {GS_Lat_E2, GD_Lat_E2}, ...
    {GS_Mean_E1, GD_Mean_E1}, {GS_Mean_E2, GD_Mean_E2}, ...
    "mean", "Pupil mean (z-score)", ...
    B_Global_Lat, B_Global_Mean, cfg);

outPngGlobal = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_Global_E1vsE2.png');
exportgraphics(figGlobal, outPngGlobal, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfGlobal = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_Global_E1vsE2.pdf');
    exportgraphics(figGlobal, outPdfGlobal, 'BackgroundColor', 'white');
end
close(figGlobal);

%% -------------------- figures: SLOPE --------------------
figLocalSlope = rbx_make_effect_figure( ...
    "Local", ...
    ["Local Standard","Local Deviant"], ...
    {LS_Lat_E1, LD_Lat_E1}, {LS_Lat_E2, LD_Lat_E2}, ...
    {LS_Slope_E1, LD_Slope_E1}, {LS_Slope_E2, LD_Slope_E2}, ...
    "Slope", char(cfg.slopeLabel), ...
    B_Local_Lat, B_Local_Slope, cfg);

outPngLocalSlope = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_LocalSlope_E1vsE2.png');
exportgraphics(figLocalSlope, outPngLocalSlope, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfLocalSlope = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_LocalSlope_E1vsE2.pdf');
    exportgraphics(figLocalSlope, outPdfLocalSlope, 'BackgroundColor', 'white');
end
close(figLocalSlope);

figGlobalSlope = rbx_make_effect_figure( ...
    "Global", ...
    ["Global Standard","Global Deviant"], ...
    {GS_Lat_E1, GD_Lat_E1}, {GS_Lat_E2, GD_Lat_E2}, ...
    {GS_Slope_E1, GD_Slope_E1}, {GS_Slope_E2, GD_Slope_E2}, ...
    "Slope", char(cfg.slopeLabel), ...
    B_Global_Lat, B_Global_Slope, cfg);

outPngGlobalSlope = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_GlobalSlope_E1vsE2.png');
exportgraphics(figGlobalSlope, outPngGlobalSlope, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfGlobalSlope = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_GlobalSlope_E1vsE2.pdf');
    exportgraphics(figGlobalSlope, outPdfGlobalSlope, 'BackgroundColor', 'white');
end
close(figGlobalSlope);

%% -------------------- outputs --------------------
OUT = struct();
OUT.subjectTable      = Tsub;
OUT.windowTable       = Twin;
OUT.statsTable        = Tstats;
OUT.windowsExperiment = Wexp;
OUT.cfg               = cfg;

OUT.files.localPNG       = outPngLocal;
OUT.files.globalPNG      = outPngGlobal;
OUT.files.localSlopePNG  = outPngLocalSlope;
OUT.files.globalSlopePNG = outPngGlobalSlope;
OUT.files.subjectCSV     = outCsv;
OUT.files.subjectMAT     = outMat;
OUT.files.windowCSV      = outWinCsv;
OUT.files.windowMAT      = outWinMat;
OUT.files.statsCSV       = outStatsCsv;
OUT.files.statsMAT       = outStatsMat;

fprintf('\n[RainBox] Saved:\n');
fprintf('  %s\n', outPngLocal);
fprintf('  %s\n', outPngGlobal);
fprintf('  %s\n', outPngLocalSlope);
fprintf('  %s\n', outPngGlobalSlope);
fprintf('  %s\n', outCsv);
fprintf('  %s\n', outWinCsv);
fprintf('  %s\n', outStatsCsv);

end

%% ========================================================================
% HELPERS
% ========================================================================

function fallbackWin = rbx_get_fallback_window(S)
if isfield(S, 'statsWinMS') && numel(S.statsWinMS)==2
    fallbackWin = double(S.statsWinMS(:)');
else
    fallbackWin = [0 2000];
end
end

function W = rbx_mask_span_window(relTimes, sigMask, fallbackWin)
relTimes = relTimes(:)';
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

function subs = rbx_list_subjects_with_file(outpath, expKey, excludeSubs)
d = dir(fullfile(outpath, 'S*'));
d = d([d.isdir]);

names = string({d.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+$', 'once'));
subs = names(keep);

if ~isempty(excludeSubs)
    subs = subs(~ismember(subs, excludeSubs));
end

hasFile = false(size(subs));
for i = 1:numel(subs)
    f = fullfile(outpath, char(subs(i)), sprintf('%s_%s_pupil_preproc.mat', subs(i), expKey));
    hasFile(i) = exist(f,'file') == 2;
end

subs = subs(hasFile);
subs = rbx_sort_subject_ids(subs);
end

function subs = rbx_sort_subject_ids(subs)
subs = string(subs);
if isempty(subs), return; end

nums = nan(size(subs));
for i = 1:numel(subs)
    tok = regexp(char(subs(i)), '^S(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        nums(i) = str2double(tok{1});
    end
end

[~,ix] = sort(nums);
subs = subs(ix);
end

function C = rbx_load_subject_curves(matFile, useField)
C = [];

S = load(matFile, 'dataTableClean', 'relTimes');
if ~isfield(S,'dataTableClean') || isempty(S.dataTableClean) || height(S.dataTableClean)==0
    return;
end
if ~isfield(S,'relTimes') || isempty(S.relTimes)
    return;
end

T = S.dataTableClean;
relTimes = S.relTimes(:)';

C = struct();
C.relTimes = relTimes;

C.LS = rbx_subject_mean_curve(T, useField, 'sound', "Local Standard");
C.LD = rbx_subject_mean_curve(T, useField, 'sound', "Local Deviant");

typeVec = string(T.type);
effectVec = strings(height(T),1);
effectVec(startsWith(typeVec, "Frequent")) = "Global Standard";
effectVec(startsWith(typeVec, "Rare"))     = "Global Deviant";

T2 = T;
T2.effect = effectVec;

C.GS = rbx_subject_mean_curve(T2, useField, 'effect', "Global Standard");
C.GD = rbx_subject_mean_curve(T2, useField, 'effect', "Global Deviant");
end

function m = rbx_subject_mean_curve(T, fieldName, condVar, condValue)
fieldName = char(fieldName);
condVec = string(T.(condVar));
idx = (condVec == string(condValue));

if ~any(idx)
    if iscell(T.(fieldName)) && ~isempty(T.(fieldName))
        m = nan(1, numel(T.(fieldName){1}));
    else
        m = nan(1,1);
    end
    return;
end

trialCells = T.(fieldName)(idx);
X = vertcat(trialCells{:});
m = mean(X, 1, 'omitnan');
end

function val = rbx_mean_in_window(relTimes, curve, startMS, endMS)
idx = relTimes >= startMS & relTimes <= endMS;
if ~any(idx) || all(isnan(curve(idx)))
    val = NaN;
else
    val = mean(curve(idx), 'omitnan');
end
end

function lat = rbx_curve_peak_latency(relTimes, curve, startMS, endMS)
idx = relTimes >= startMS & relTimes <= endMS;
if ~any(idx)
    lat = NaN;
    return;
end

t = relTimes(idx);
y = curve(idx);

good = isfinite(t) & isfinite(y);
t = t(good);
y = y(good);

if isempty(y)
    lat = NaN;
    return;
end

[~,ix] = max(y);
lat = t(ix);
end

function slopeVal = rbx_curve_slope(relTimes, curve, startMS, endMS, slopeMode, slopeUnit)
% Default/recommended:
%   slope from window START to subject-specific PEAK within that window

idx = relTimes >= startMS & relTimes <= endMS;
if ~any(idx)
    slopeVal = NaN;
    return;
end

t = relTimes(idx);
y = curve(idx);

good = isfinite(t) & isfinite(y);
t = t(good);
y = y(good);

if numel(t) < 2
    slopeVal = NaN;
    return;
end

switch lower(char(slopeMode))
    case 'full_window'
        tFit = t;
        yFit = y;

    otherwise % 'start_to_peak'
        [~,ixPeak] = max(y);
        tFit = t(1:ixPeak);
        yFit = y(1:ixPeak);
end

if numel(tFit) < 2 || numel(unique(tFit)) < 2
    slopeVal = NaN;
    return;
end

switch lower(char(slopeUnit))
    case 'ms'
        x = tFit - tFit(1);         % ms
    otherwise
        x = (tFit - tFit(1)) / 1000; % seconds
end

p = polyfit(x, yFit, 1);
slopeVal = p(1);
end

function Twin = rbx_build_window_table(Wexp, cfg)
Effect = [ ...
    "Local"; "Local"; ...
    "Global"; "Global"];

Experiment = [ ...
    cfg.labelE1; cfg.labelE2; ...
    cfg.labelE1; cfg.labelE2];

StartMS = [ ...
    Wexp.Local.E1.startMS; ...
    Wexp.Local.E2.startMS; ...
    Wexp.Global.E1.startMS; ...
    Wexp.Global.E2.startMS];

EndMS = [ ...
    Wexp.Local.E1.endMS; ...
    Wexp.Local.E2.endMS; ...
    Wexp.Global.E1.endMS; ...
    Wexp.Global.E2.endMS];

UsedFallback = [ ...
    Wexp.Local.E1.usedFallback; ...
    Wexp.Local.E2.usedFallback; ...
    Wexp.Global.E1.usedFallback; ...
    Wexp.Global.E2.usedFallback];

UsedForConditions = [ ...
    "Local Standard; Local Deviant"; ...
    "Local Standard; Local Deviant"; ...
    "Global Standard; Global Deviant"; ...
    "Global Standard; Global Deviant"];

WindowType = repmat("Experiment-specific blue bar", 4, 1);

Twin = table(Effect, WindowType, Experiment, StartMS, EndMS, UsedFallback, UsedForConditions);
end

function Tstats = rbx_build_stats_table(cfg, ...
    LS_Lat_E1, LS_Lat_E2, LD_Lat_E1, LD_Lat_E2, ...
    LS_Mean_E1, LS_Mean_E2, LD_Mean_E1, LD_Mean_E2, ...
    LS_Slope_E1, LS_Slope_E2, LD_Slope_E1, LD_Slope_E2, ...
    GS_Lat_E1, GS_Lat_E2, GD_Lat_E1, GD_Lat_E2, ...
    GS_Mean_E1, GS_Mean_E2, GD_Mean_E1, GD_Mean_E2, ...
    GS_Slope_E1, GS_Slope_E2, GD_Slope_E1, GD_Slope_E2)

Effect          = strings(18,1);
Measure         = strings(18,1);
Condition       = strings(18,1);
Family          = strings(18,1);
ExperimentA     = repmat(cfg.labelE1, 18, 1);
ExperimentB     = repmat(cfg.labelE2, 18, 1);
Npairs          = nan(18,1);
RawPValue       = nan(18,1);
CorrectedPValue = nan(18,1);
Significance    = strings(18,1);

row = 0;

% ---------- LOCAL LATENCY ----------
row=row+1; [p,~,n] = rbx_paired_pvalue(LS_Lat_E1, LS_Lat_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Peak latency"; Condition(row)="Local Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(LD_Lat_E1, LD_Lat_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Peak latency"; Condition(row)="Local Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(LD_Lat_E1 - LS_Lat_E1, LD_Lat_E2 - LS_Lat_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Peak latency"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- LOCAL MEAN ----------
row=row+1; [p,~,n] = rbx_paired_pvalue(LS_Mean_E1, LS_Mean_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Mean pupil size"; Condition(row)="Local Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(LD_Mean_E1, LD_Mean_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Mean pupil size"; Condition(row)="Local Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(LD_Mean_E1 - LS_Mean_E1, LD_Mean_E2 - LS_Mean_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Mean pupil size"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- LOCAL SLOPE ----------
row=row+1; [p,~,n] = rbx_paired_pvalue(LS_Slope_E1, LS_Slope_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Slope pupil size"; Condition(row)="Local Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(LD_Slope_E1, LD_Slope_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Slope pupil size"; Condition(row)="Local Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(LD_Slope_E1 - LS_Slope_E1, LD_Slope_E2 - LS_Slope_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Slope pupil size"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;


% ---------- GLOBAL LATENCY ----------
row=row+1; [p,~,n] = rbx_paired_pvalue(GS_Lat_E1, GS_Lat_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Peak latency"; Condition(row)="Global Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(GD_Lat_E1, GD_Lat_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Peak latency"; Condition(row)="Global Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(GD_Lat_E1 - GS_Lat_E1, GD_Lat_E2 - GS_Lat_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Peak latency"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- GLOBAL MEAN ----------
row=row+1; [p,~,n] = rbx_paired_pvalue(GS_Mean_E1, GS_Mean_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Mean pupil size"; Condition(row)="Global Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(GD_Mean_E1, GD_Mean_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Mean pupil size"; Condition(row)="Global Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(GD_Mean_E1 - GS_Mean_E1, GD_Mean_E2 - GS_Mean_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Mean pupil size"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- GLOBAL SLOPE ----------
row=row+1; [p,~,n] = rbx_paired_pvalue(GS_Slope_E1, GS_Slope_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Slope pupil size"; Condition(row)="Global Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(GD_Slope_E1, GD_Slope_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Slope pupil size"; Condition(row)="Global Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = rbx_paired_pvalue(GD_Slope_E1 - GS_Slope_E1, GD_Slope_E2 - GS_Slope_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Slope pupil size"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- define family labels ----------
for k = 1:numel(Effect)
    Family(k) = Effect(k) + " | " + Measure(k);
end

% ---------- correction WITHIN each 3-test family ----------
uFam = unique(Family, 'stable');
for iFam = 1:numel(uFam)
    idxFam = (Family == uFam(iFam));
    CorrectedPValue(idxFam) = rbx_correct_pvalues(RawPValue(idxFam), cfg.multCompMethod);
end

% ---------- stars from corrected p-values ----------
for k = 1:numel(CorrectedPValue)
    Significance(k) = rbx_p_to_stars(CorrectedPValue(k));
end

Tstats = table( ...
    Effect, Measure, Family, Condition, ExperimentA, ExperimentB, ...
    Npairs, RawPValue, CorrectedPValue, Significance);
end

function fig = rbx_make_effect_figure(effectName, condLabels, latE1Cell, latE2Cell, rightE1Cell, rightE2Cell, rightPanelTitle, rightYLabel, BLat, BRight, cfg)

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 cfg.figW cfg.figH], ...
    'InvertHardcopy','off');

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
tl.OuterPosition = [0.02 0.04 0.98 0.90];

groupCenters = [1 2];
xA = groupCenters - 0.15;
xB = groupCenters + 0.15;

%% ---------- LEFT: LATENCY ----------
ax1 = nexttile(tl); hold(ax1,'on');
rbx_style_axis(ax1, cfg.bgLatency);

title(ax1, sprintf('%s effects | Pupil peak latency', effectName), 'FontSize', 12, 'FontWeight','bold');
ylabel(ax1, 'Pupil latency (ms)', 'FontSize', 12);

latSupport = [];
latDataAll = [];

for c = 1:2
    y1 = latE1Cell{c};
    y2 = latE2Cell{c};

    [v1min,v1max] = rbx_draw_half_violin(ax1, y1, xA(c), 'left',  cfg.colorE1, cfg.rainWidth, cfg.violinAlpha);
    [v2min,v2max] = rbx_draw_half_violin(ax1, y2, xB(c), 'right', cfg.colorE2, cfg.rainWidth, cfg.violinAlpha);

    latSupport = [latSupport; v1min; v1max; v2min; v2max]; %#ok<AGROW>
    latDataAll = [latDataAll; y1(:); y2(:)]; %#ok<AGROW>

    rbx_draw_paired_lines(ax1, y1, y2, xA(c), xB(c), cfg.grayLine);
    rbx_draw_custom_box(ax1, y1, xA(c), cfg.colorE1, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    rbx_draw_custom_box(ax1, y2, xB(c), cfg.colorE2, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    rbx_draw_points(ax1, y1, xA(c), cfg.colorE1, cfg.jitter, cfg.pointSize);
    rbx_draw_points(ax1, y2, xB(c), cfg.colorE2, cfg.jitter, cfg.pointSize);
end

set(ax1, 'XLim', [0.5 2.5], 'XTick', groupCenters, 'XTickLabel', cellstr(condLabels));
rbx_set_nice_ylim(ax1, latDataAll, latSupport);

rbx_add_pair_bracket(ax1, xA(1), xB(1), [latE1Cell{1}; latE2Cell{1}], 1, BLat.cond1);
rbx_add_pair_bracket(ax1, xA(2), xB(2), [latE1Cell{2}; latE2Cell{2}], 2, BLat.cond2);
rbx_add_global_bracket(ax1, groupCenters(1), groupCenters(2), latDataAll, BLat.contrast);

%% ---------- RIGHT: GENERIC PARAMETER ----------
ax2 = nexttile(tl); hold(ax2,'on');
rbx_style_axis(ax2, cfg.bgPupil);

title(ax2, sprintf('%s effects | Pupil %s size', effectName, lower(rightPanelTitle)), 'FontSize', 12, 'FontWeight','bold');
ylabel(ax2, rightYLabel, 'FontSize', 12);

rightSupport = [];
rightDataAll = [];

for c = 1:2
    y1 = rightE1Cell{c};
    y2 = rightE2Cell{c};

    [v1min,v1max] = rbx_draw_half_violin(ax2, y1, xA(c), 'left',  cfg.colorE1, cfg.rainWidth, cfg.violinAlpha);
    [v2min,v2max] = rbx_draw_half_violin(ax2, y2, xB(c), 'right', cfg.colorE2, cfg.rainWidth, cfg.violinAlpha);

    rightSupport = [rightSupport; v1min; v1max; v2min; v2max]; %#ok<AGROW>
    rightDataAll = [rightDataAll; y1(:); y2(:)]; %#ok<AGROW>

    rbx_draw_paired_lines(ax2, y1, y2, xA(c), xB(c), cfg.grayLine);
    rbx_draw_custom_box(ax2, y1, xA(c), cfg.colorE1, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    rbx_draw_custom_box(ax2, y2, xB(c), cfg.colorE2, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    rbx_draw_points(ax2, y1, xA(c), cfg.colorE1, cfg.jitter, cfg.pointSize);
    rbx_draw_points(ax2, y2, xB(c), cfg.colorE2, cfg.jitter, cfg.pointSize);
end

set(ax2, 'XLim', [0.5 2.5], 'XTick', groupCenters, 'XTickLabel', cellstr(condLabels));
rbx_set_nice_ylim(ax2, rightDataAll, rightSupport);

rbx_add_pair_bracket(ax2, xA(1), xB(1), [rightE1Cell{1}; rightE2Cell{1}], 1, BRight.cond1);
rbx_add_pair_bracket(ax2, xA(2), xB(2), [rightE1Cell{2}; rightE2Cell{2}], 2, BRight.cond2);
rbx_add_global_bracket(ax2, groupCenters(1), groupCenters(2), rightDataAll, BRight.contrast);

%% ---------- Legend in LEFT plot ----------
h1L = plot(ax1, nan, nan, '-', 'Color', cfg.colorE1, 'LineWidth', 3, 'DisplayName', char(cfg.labelE1));
h2L = plot(ax1, nan, nan, '-', 'Color', cfg.colorE2, 'LineWidth', 3, 'DisplayName', char(cfg.labelE2));

lgd1 = legend(ax1, [h1L h2L], 'Location', 'none', 'Box', 'off');
lgd1.FontSize = 10;
lgd1.Units = 'normalized';

ax1Pos = ax1.Position;
legendW = 0.22;
legendH = 0.08;
legendX = ax1Pos(1) + ax1Pos(3)/2 - legendW/2;
legendY = ax1Pos(2) + 0.03;
lgd1.Position = [legendX legendY legendW legendH];

%% ---------- Legend in RIGHT plot ----------
h1R = plot(ax2, nan, nan, '-', 'Color', cfg.colorE1, 'LineWidth', 3, 'DisplayName', char(cfg.labelE1));
h2R = plot(ax2, nan, nan, '-', 'Color', cfg.colorE2, 'LineWidth', 3, 'DisplayName', char(cfg.labelE2));

lgd2 = legend(ax2, [h1R h2R], 'Location', 'none', 'Box', 'off');
lgd2.FontSize = 10;
lgd2.Units = 'normalized';

ax2Pos = ax2.Position;
legendX = ax2Pos(1) + ax2Pos(3)/2 - legendW/2;
legendY = ax2Pos(2) + 0.03;
lgd2.Position = [legendX legendY legendW legendH];

sgtitle(fig, sprintf('Pupil size | %s effects', effectName), 'FontSize', 12, 'FontWeight', 'bold');
end

function rbx_style_axis(ax, bgColor)
set(ax, ...
    'Color', bgColor, ...
    'FontSize', 11, ...
    'LineWidth', 1, ...
    'Box', 'off', ...
    'Layer', 'top', ...
    'TickDir', 'out');
grid(ax, 'off');
end

function [vmin, vmax] = rbx_draw_half_violin(ax, y, x0, sideStr, rgb, maxWidth, faceAlpha)
y = y(isfinite(y));

if numel(y) < 2
    vmin = nan;
    vmax = nan;
    return;
end

try
    if numel(unique(y)) == 1
        yy = linspace(y(1)-eps(y(1)+1), y(1)+eps(y(1)+1), 64);
        f  = ones(size(yy));
    else
        [f, yy] = ksdensity(y, 'NumPoints', 200);
    end
catch
    [N, edges] = histcounts(y, 30, 'Normalization', 'pdf');
    yy = edges(1:end-1) + diff(edges)/2;
    f = N;
end

if isempty(f) || all(~isfinite(f)) || max(f) <= 0
    vmin = nan;
    vmax = nan;
    return;
end

f = f ./ max(f);
f = f * maxWidth;

switch lower(char(sideStr))
    case 'left'
        xCurve = x0 - f;
    otherwise
        xCurve = x0 + f;
end

xPoly = [repmat(x0,1,numel(yy)) xCurve(end:-1:1)];
yPoly = [yy yy(end:-1:1)];

patch(ax, xPoly, yPoly, rgb, ...
    'FaceAlpha', faceAlpha, ...
    'EdgeColor', rgb, ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

vmin = min(yy);
vmax = max(yy);
end

function rbx_draw_paired_lines(ax, y1, y2, x1, x2, lineColor)
m = isfinite(y1) & isfinite(y2);
idx = find(m);
for i = 1:numel(idx)
    ii = idx(i);
    plot(ax, [x1 x2], [y1(ii) y2(ii)], '-', ...
        'Color', lineColor, 'LineWidth', 1.5, 'HandleVisibility','off');
end
end

function rbx_draw_points(ax, y, x0, rgb, jitterAmt, markerSize)
y = y(:);
y = y(isfinite(y));
if isempty(y), return; end

rng(1);
x = x0 + (rand(size(y)) - 0.5) * 2 * jitterAmt;

scatter(ax, x, y, markerSize, ...
    'MarkerFaceColor', rgb, ...
    'MarkerEdgeColor', [0.30 0.30 0.30], ...
    'LineWidth', 0.35, ...
    'MarkerFaceAlpha', 0.35, ...
    'MarkerEdgeAlpha', 0.35, ...
    'HandleVisibility', 'off');
end

function rbx_draw_custom_box(ax, y, x0, rgb, boxW, boxLW, meanLW)
y = y(isfinite(y));
if isempty(y), return; end

q  = prctile(y, [25 75]);
q1 = q(1);
q3 = q(2);
mu = mean(y, 'omitnan');

iqrVal  = q3 - q1;
loFence = q1 - 1.5 * iqrVal;
hiFence = q3 + 1.5 * iqrVal;

yIn = y(y >= loFence & y <= hiFence);
if isempty(yIn)
    yLow  = min(y);
    yHigh = max(y);
else
    yLow  = min(yIn);
    yHigh = max(yIn);
end

rectangle(ax, 'Position', [x0 - boxW/2, q1, boxW, q3 - q1], ...
    'EdgeColor', rgb, 'LineWidth', boxLW, 'FaceColor', 'none');

plot(ax, [x0 x0], [yLow q1],  '-', 'Color', rgb, 'LineWidth', boxLW, 'HandleVisibility','off');
plot(ax, [x0 x0], [q3 yHigh], '-', 'Color', rgb, 'LineWidth', boxLW, 'HandleVisibility','off');

capW = boxW * 0.45;
plot(ax, [x0-capW/2 x0+capW/2], [yLow yLow],   '-', 'Color', rgb, 'LineWidth', boxLW, 'HandleVisibility','off');
plot(ax, [x0-capW/2 x0+capW/2], [yHigh yHigh], '-', 'Color', rgb, 'LineWidth', boxLW, 'HandleVisibility','off');

plot(ax, [x0-boxW/2 x0+boxW/2], [mu mu], '-', ...
    'Color', rgb, 'LineWidth', meanLW, 'HandleVisibility','off');
end

function rbx_set_nice_ylim(ax, yData, ySupport)
yData = yData(isfinite(yData));
ySupport = ySupport(isfinite(ySupport));

if isempty(yData) && isempty(ySupport)
    return;
end

if isempty(ySupport)
    yMin = min(yData);
    yMax = max(yData);
else
    yMin = min([yData; ySupport]);
    yMax = max([yData; ySupport]);
end

yr = yMax - yMin;
if yr <= 0
    yr = max(abs(yMax), 1);
end

padLo = 0.10 * yr;
padHi = 0.15 * yr;

ylim(ax, [yMin - padLo, yMax + padHi]);
end

function yTop = rbx_add_pair_bracket(ax, x1, x2, yVals, groupIdx, txt)
yl = ylim(ax);
yr = yl(2) - yl(1);

yVals = yVals(isfinite(yVals));
if isempty(yVals)
    yTop = NaN;
    return;
end

pairTop = max(yVals);

if groupIdx == 1
    y = min(pairTop + 0.08*yr, yl(2) - 0.18*yr);
else
    y = min(pairTop + 0.06*yr, yl(2) - 0.28*yr);
end

h = 0.025 * yr;
txtGap = 0.015 * yr;

plot(ax, [x1 x1 x2 x2], [y y+h y+h y], '-', ...
    'Color', 'k', 'LineWidth', 1.6, 'HandleVisibility','off');

text(ax, mean([x1 x2]), y+h+txtGap, char(txt), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'Color', 'k', ...
    'HandleVisibility', 'off');

yTop = y + h + txtGap;
end

function rbx_add_global_bracket(ax, x1, x2, ~, txt)
yl = ylim(ax);
yr = yl(2) - yl(1);

% Fixed height relative to current axis top
y = yl(2) - 0.12 * yr;
h = 0.020 * yr;
txtGap = 0.008 * yr;

plot(ax, [x1 x1 x2 x2], [y y+h y+h y], '-', ...
    'Color', 'k', 'LineWidth', 1.6, 'HandleVisibility','off');

text(ax, mean([x1 x2]), y+h+txtGap, char(txt), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'Color', 'k', ...
    'HandleVisibility', 'off');
end

function [p, stars, n] = rbx_paired_pvalue(a, b, methodStr)
m = isfinite(a) & isfinite(b);
a = a(m);
b = b(m);
n = numel(a);

if n < 2
    p = NaN;
    stars = "n/a";
    return;
end

if all(abs(a-b) < 1e-12)
    p = 1;
    stars = "ns";
    return;
end

switch lower(char(methodStr))
    case 'ttest'
        [~,p] = ttest(a,b);
    otherwise
        try
            p = signrank(a,b);
        catch
            [~,p] = ttest(a,b);
        end
end

stars = rbx_p_to_stars(p);
end

function txt = rbx_p_to_stars(p)
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

function B = rbx_get_bracket_info(Tstats, effectName, measureName)
idx = (Tstats.Effect == string(effectName)) & (Tstats.Measure == string(measureName));
Ts = Tstats(idx,:);

if height(Ts) ~= 3
    error('Expected exactly 3 rows for %s | %s in Tstats.', effectName, measureName);
end

B = struct();
B.cond1    = Ts.Significance(1);
B.cond2    = Ts.Significance(2);
B.contrast = Ts.Significance(3);

B.p1 = Ts.CorrectedPValue(1);
B.p2 = Ts.CorrectedPValue(2);
B.p3 = Ts.CorrectedPValue(3);
end

function pAdj = rbx_correct_pvalues(pRaw, methodStr)
pAdj = nan(size(pRaw));

ok = isfinite(pRaw);
if ~any(ok)
    return;
end

p = pRaw(ok);

switch lower(char(methodStr))
    case 'none'
        pAdj(ok) = p;

    case 'bonferroni'
        pAdj(ok) = min(p * numel(p), 1);

    otherwise % 'holm'
        pAdj(ok) = rbx_holm_adjust(p);
end
end

function pAdj = rbx_holm_adjust(p)
origSize = size(p);
p = p(:);

m = numel(p);
[ps, ord] = sort(p, 'ascend');

adjSorted = (m - (1:m)' + 1) .* ps;
adjSorted = min(adjSorted, 1);
adjSorted = cummax(adjSorted);   % enforce monotonicity

pAdj = nan(size(p));
pAdj(ord) = adjSorted;
pAdj = reshape(pAdj, origSize);
end