function OUT = eeg_export_rainbox_compare_E1E2_Pz(outroot, grandFold, GROUPPMAP, excludeSubs, targets, chanName)
% eeg_export_rainbox_compare_E1E2_Pz
%
% EEG rainbox summary plots for E1 vs E2 at ONE channel (default: Pz),
% mirroring the logic/style of your pupil rainbox function.
%
% Uses:
%   - auditory_all_LG_Sxx_E1.set
%   - auditory_all_LG_Sxx_E2.set
%
% Reconstructs 4 condition means per subject from eventlg_type:
%   FLS = Frequent Local Standard
%   RLS = Rare Local Standard
%   FLD = Frequent Local Deviant
%   RLD = Rare Local Deviant
%
% Then derives (equal-weight factorial averages):
%   Local Standard  = mean(FLS, RLS)
%   Local Deviant   = mean(FLD, RLD)
%   Global Standard = mean(FLS, FLD)
%   Global Deviant  = mean(RLS, RLD)
%
% Summary metrics per subject:
%   - peak latency (ms)
%   - mean amplitude (uV)
%   - slope from window start to subject-specific peak (uV/s by default)
%
% Window logic:
%   - Uses the FULL blue-bar span from GROUPPMAP.(E1/E2).mLoc and .mGlo
%   - Same experiment-specific window for Standard and Deviant of that effect
%
% Multiple-comparison correction:
%   - Holm within each Effect x Measure 3-test family:
%       Standard, Deviant, Deviant-Standard contrast
%
% Outputs:
%   1) EEG_RAINBOX_COMPARE_Pz_Local_E1vsE2.png
%   2) EEG_RAINBOX_COMPARE_Pz_Global_E1vsE2.png
%   3) EEG_RAINBOX_COMPARE_Pz_LocalSlope_E1vsE2.png
%   4) EEG_RAINBOX_COMPARE_Pz_GlobalSlope_E1vsE2.png
%   + subject/window/stats tables
%
% Miguel Contreras Altamirano, 2026
% Adapted for EEG Pz rainbox export

%% -------------------- defaults --------------------
if nargin < 4 || isempty(excludeSubs)
    excludeSubs = strings(0,1);
end
if nargin < 5 || isempty(targets)
    targets = [ ...
        "Frequent Local Standard", ...
        "Rare Local Standard", ...
        "Frequent Local Deviant", ...
        "Rare Local Deviant"];
end
if nargin < 6 || isempty(chanName)
    chanName = "Pz";
end

excludeSubs = string(excludeSubs);
targets     = string(targets);
chanName    = string(chanName);

if ~exist(outroot, 'dir')
    error('outroot does not exist: %s', outroot);
end
if ~exist(grandFold, 'dir')
    mkdir(grandFold);
end

if exist('pop_loadset','file') ~= 2
    error('pop_loadset not found on MATLAB path. Make sure EEGLAB is loaded first.');
end

%% -------------------- style/config --------------------
cfg = struct();

cfg.labelE1 = "Counting";
cfg.labelE2 = "Visual interference";
cfg.chanName = chanName;

cfg.colorE1 = [0 0.45 0.85];
cfg.colorE2 = [0.90 0.45 0.10];
cfg.grayLine = [0.62 0.62 0.62];

% cfg.bgLatency = [0.93 0.93 0.93]; % background color left panel
% cfg.bgERP     = [0.99 0.94 0.89]; % background color right panel

cfg.bgLatency = [1 1 1];
cfg.bgERP     = [1 1 1];

cfg.statsMethod = "signrank";   % "signrank" | "ttest"
cfg.multCompMethod = "holm";    % "holm" | "bonferroni" | "none"

cfg.dpi     = 300;
cfg.savePDF = false;

cfg.figW = 13.8;
cfg.figH = 4.8;

cfg.rainWidth     = 0.13;
cfg.boxWidth      = 0.18;
cfg.jitter        = 0.022;
cfg.pointSize     = 28;
cfg.boxLineWidth  = 2.0;
cfg.meanLineWidth = 3.0;
cfg.violinAlpha   = 0.22;

% Peak/slope settings
cfg.peakPolarity = "positive";   % <---- change to "negative" if ever needed
cfg.slopeMode    = "start_to_peak";  % "start_to_peak" | "full_window"
cfg.slopeUnit    = "ms";              % "s" or "ms"

% Labels title
cfg.measureLatency = "Peak latency";
cfg.measureMean    = "Mean amplitude";
cfg.measureSlope   = "Slope amplitude";

fprintf(['\n[EEG RainBox] CHANNEL = %s\n'], char(cfg.chanName));
fprintf(['[EEG RainBox] MULTIPLE-COMPARISON CORRECTION = %s within each ' ...
    'Effect x Measure family (3 tests: Standard, Deviant, Deviant-Standard contrast)\n'], ...
    upper(char(cfg.multCompMethod)));
fprintf('[EEG RainBox] PEAK POLARITY = %s\n', upper(char(cfg.peakPolarity)));
fprintf('[EEG RainBox] SLOPE DEFINITION = linear slope from window start to subject peak inside that window (%s)\n', cfg.slopeUnit);

%% -------------------- sanity checks --------------------
mustHaveFields = {'E1','E2'};
for k = 1:numel(mustHaveFields)
    if ~isfield(GROUPPMAP, mustHaveFields{k})
        error('GROUPPMAP.%s is missing.', mustHaveFields{k});
    end
end

req = {'times','mLoc','mGlo'};
for k = 1:numel(req)
    if ~isfield(GROUPPMAP.E1, req{k})
        error('GROUPPMAP.E1.%s is missing.', req{k});
    end
    if ~isfield(GROUPPMAP.E2, req{k})
        error('GROUPPMAP.E2.%s is missing.', req{k});
    end
end

%% -------------------- windows from EEG blue bars --------------------

fallbackE1 = erbx_get_fallback_window(GROUPPMAP.E1.times);
fallbackE2 = erbx_get_fallback_window(GROUPPMAP.E2.times);

Wexp.Local.E1  = erbx_mask_span_window(GROUPPMAP.E1.times, GROUPPMAP.E1.mLoc, fallbackE1);
Wexp.Local.E2  = erbx_mask_span_window(GROUPPMAP.E2.times, GROUPPMAP.E2.mLoc, fallbackE2);

Wexp.Global.E1 = erbx_mask_span_window(GROUPPMAP.E1.times, GROUPPMAP.E1.mGlo, fallbackE1);
Wexp.Global.E2 = erbx_mask_span_window(GROUPPMAP.E2.times, GROUPPMAP.E2.mGlo, fallbackE2);

% ---------------------------------------------------------------------
% NEW: if E2 has no significant mask, reuse the corresponding E1 window
% ---------------------------------------------------------------------
if Wexp.Local.E2.usedFallback
    Wexp.Local.E2.startMS = Wexp.Local.E1.startMS;
    Wexp.Local.E2.endMS   = Wexp.Local.E1.endMS;
    Wexp.Local.E2.idx     = GROUPPMAP.E2.times >= Wexp.Local.E2.startMS & ...
                            GROUPPMAP.E2.times <= Wexp.Local.E2.endMS;
    Wexp.Local.E2.borrowedFromE1 = true;
else
    Wexp.Local.E2.borrowedFromE1 = false;
end

if Wexp.Global.E2.usedFallback
    Wexp.Global.E2.startMS = Wexp.Global.E1.startMS;
    Wexp.Global.E2.endMS   = Wexp.Global.E1.endMS;
    Wexp.Global.E2.idx     = GROUPPMAP.E2.times >= Wexp.Global.E2.startMS & ...
                             GROUPPMAP.E2.times <= Wexp.Global.E2.endMS;
    Wexp.Global.E2.borrowedFromE1 = true;
else
    Wexp.Global.E2.borrowedFromE1 = false;
end


fprintf('\n[EEG RainBox] WINDOWS USED FOR EXTRACTION\n');
fprintf('  Local  | %s: %.0f..%.0f ms | %s: %.0f..%.0f ms\n', ...
    cfg.labelE1, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, ...
    cfg.labelE2, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);

fprintf('  Global | %s: %.0f..%.0f ms | %s: %.0f..%.0f ms\n', ...
    cfg.labelE1, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, ...
    cfg.labelE2, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);

fprintf('  Local  | E2 borrowed E1 window: %d\n', Wexp.Local.E2.borrowedFromE1);
fprintf('  Global | E2 borrowed E1 window: %d\n', Wexp.Global.E2.borrowedFromE1);


%% -------------------- subjects present in BOTH E1 and E2 --------------------
subsE1 = erbx_list_subjects_with_file(outroot, 'E1', excludeSubs);
subsE2 = erbx_list_subjects_with_file(outroot, 'E2', excludeSubs);

subNames = subsE1(ismember(subsE1, subsE2));
subNames = erbx_sort_subject_ids(subNames);

if isempty(subNames)
    error('No common subjects found with both E1 and E2 EEG .set files.');
end

%% -------------------- subject-level extraction --------------------
nSub = numel(subNames);

LS_Lat_E1 = nan(nSub,1); LS_Lat_E2 = nan(nSub,1);
LD_Lat_E1 = nan(nSub,1); LD_Lat_E2 = nan(nSub,1);

LS_Mean_E1 = nan(nSub,1); LS_Mean_E2 = nan(nSub,1);
LD_Mean_E1 = nan(nSub,1); LD_Mean_E2 = nan(nSub,1);

LS_Slope_E1 = nan(nSub,1); LS_Slope_E2 = nan(nSub,1);
LD_Slope_E1 = nan(nSub,1); LD_Slope_E2 = nan(nSub,1);

GS_Lat_E1 = nan(nSub,1); GS_Lat_E2 = nan(nSub,1);
GD_Lat_E1 = nan(nSub,1); GD_Lat_E2 = nan(nSub,1);

GS_Mean_E1 = nan(nSub,1); GS_Mean_E2 = nan(nSub,1);
GD_Mean_E1 = nan(nSub,1); GD_Mean_E2 = nan(nSub,1);

GS_Slope_E1 = nan(nSub,1); GS_Slope_E2 = nan(nSub,1);
GD_Slope_E1 = nan(nSub,1); GD_Slope_E2 = nan(nSub,1);

for i = 1:nSub
    sID = subNames(i);

    fileE1 = fullfile(outroot, char(sID), sprintf('auditory_all_LG_%s_E1.set', char(sID)));
    fileE2 = fullfile(outroot, char(sID), sprintf('auditory_all_LG_%s_E2.set', char(sID)));

    if ~exist(fileE1,'file') || ~exist(fileE2,'file')
        continue;
    end

    C1 = erbx_load_subject_curves_eeg(fileE1, cfg.chanName, targets);
    C2 = erbx_load_subject_curves_eeg(fileE2, cfg.chanName, targets);

    if isempty(C1) || isempty(C2)
        continue;
    end

    % ---------- LOCAL ----------
    LS_Mean_E1(i) = erbx_mean_in_window(C1.relTimes, C1.LS, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS);
    LD_Mean_E1(i) = erbx_mean_in_window(C1.relTimes, C1.LD, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS);

    LS_Mean_E2(i) = erbx_mean_in_window(C2.relTimes, C2.LS, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);
    LD_Mean_E2(i) = erbx_mean_in_window(C2.relTimes, C2.LD, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS);

    LS_Lat_E1(i) = erbx_curve_peak_latency(C1.relTimes, C1.LS, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, cfg.peakPolarity);
    LD_Lat_E1(i) = erbx_curve_peak_latency(C1.relTimes, C1.LD, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, cfg.peakPolarity);

    LS_Lat_E2(i) = erbx_curve_peak_latency(C2.relTimes, C2.LS, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS, cfg.peakPolarity);
    LD_Lat_E2(i) = erbx_curve_peak_latency(C2.relTimes, C2.LD, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS, cfg.peakPolarity);

    LS_Slope_E1(i) = erbx_curve_slope(C1.relTimes, C1.LS, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);
    LD_Slope_E1(i) = erbx_curve_slope(C1.relTimes, C1.LD, Wexp.Local.E1.startMS, Wexp.Local.E1.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);

    LS_Slope_E2(i) = erbx_curve_slope(C2.relTimes, C2.LS, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);
    LD_Slope_E2(i) = erbx_curve_slope(C2.relTimes, C2.LD, Wexp.Local.E2.startMS, Wexp.Local.E2.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);

    % ---------- GLOBAL ----------
    GS_Mean_E1(i) = erbx_mean_in_window(C1.relTimes, C1.GS, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);
    GD_Mean_E1(i) = erbx_mean_in_window(C1.relTimes, C1.GD, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS);

    GS_Mean_E2(i) = erbx_mean_in_window(C2.relTimes, C2.GS, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);
    GD_Mean_E2(i) = erbx_mean_in_window(C2.relTimes, C2.GD, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS);

    GS_Lat_E1(i) = erbx_curve_peak_latency(C1.relTimes, C1.GS, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, cfg.peakPolarity);
    GD_Lat_E1(i) = erbx_curve_peak_latency(C1.relTimes, C1.GD, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, cfg.peakPolarity);

    GS_Lat_E2(i) = erbx_curve_peak_latency(C2.relTimes, C2.GS, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS, cfg.peakPolarity);
    GD_Lat_E2(i) = erbx_curve_peak_latency(C2.relTimes, C2.GD, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS, cfg.peakPolarity);

    GS_Slope_E1(i) = erbx_curve_slope(C1.relTimes, C1.GS, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);
    GD_Slope_E1(i) = erbx_curve_slope(C1.relTimes, C1.GD, Wexp.Global.E1.startMS, Wexp.Global.E1.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);

    GS_Slope_E2(i) = erbx_curve_slope(C2.relTimes, C2.GS, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);
    GD_Slope_E2(i) = erbx_curve_slope(C2.relTimes, C2.GD, Wexp.Global.E2.startMS, Wexp.Global.E2.endMS, cfg.slopeMode, cfg.slopeUnit, cfg.peakPolarity);
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
    'LocalStandard_Mean_E1_uV',    'LocalStandard_Mean_E2_uV', ...
    'LocalDeviant_Mean_E1_uV',     'LocalDeviant_Mean_E2_uV', ...
    'LocalStandard_Slope_E1',      'LocalStandard_Slope_E2', ...
    'LocalDeviant_Slope_E1',       'LocalDeviant_Slope_E2', ...
    'GlobalStandard_Latency_E1_ms','GlobalStandard_Latency_E2_ms', ...
    'GlobalDeviant_Latency_E1_ms', 'GlobalDeviant_Latency_E2_ms', ...
    'GlobalStandard_Mean_E1_uV',   'GlobalStandard_Mean_E2_uV', ...
    'GlobalDeviant_Mean_E1_uV',    'GlobalDeviant_Mean_E2_uV', ...
    'GlobalStandard_Slope_E1',     'GlobalStandard_Slope_E2', ...
    'GlobalDeviant_Slope_E1',      'GlobalDeviant_Slope_E2'});

outCsv = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_subject_values.csv', char(cfg.chanName)));
outMat = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_subject_values.mat', char(cfg.chanName)));
writetable(Tsub, outCsv);
save(outMat, 'Tsub', 'Wexp', 'cfg', '-v7.3');

%% -------------------- blue-bar windows table --------------------
Twin = erbx_build_window_table(Wexp, cfg);

outWinCsv = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_sig_windows.csv', char(cfg.chanName)));
outWinMat = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_sig_windows.mat', char(cfg.chanName)));
writetable(Twin, outWinCsv);
save(outWinMat, 'Twin', 'Wexp', '-v7.3');

disp(Twin);

%% -------------------- stats summary table --------------------
Tstats = erbx_build_stats_table(cfg, ...
    LS_Lat_E1, LS_Lat_E2, LD_Lat_E1, LD_Lat_E2, ...
    LS_Mean_E1, LS_Mean_E2, LD_Mean_E1, LD_Mean_E2, ...
    LS_Slope_E1, LS_Slope_E2, LD_Slope_E1, LD_Slope_E2, ...
    GS_Lat_E1, GS_Lat_E2, GD_Lat_E1, GD_Lat_E2, ...
    GS_Mean_E1, GS_Mean_E2, GD_Mean_E1, GD_Mean_E2, ...
    GS_Slope_E1, GS_Slope_E2, GD_Slope_E1, GD_Slope_E2);

outStatsCsv = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_stats_summary.csv', char(cfg.chanName)));
outStatsMat = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_stats_summary.mat', char(cfg.chanName)));
writetable(Tstats, outStatsCsv);
save(outStatsMat, 'Tstats', '-v7.3');

disp(Tstats);

%% -------------------- bracket labels --------------------
B_Local_Lat   = erbx_get_bracket_info(Tstats, "Local",  cfg.measureLatency);
B_Local_Mean  = erbx_get_bracket_info(Tstats, "Local",  cfg.measureMean);
B_Local_Slope = erbx_get_bracket_info(Tstats, "Local",  cfg.measureSlope);

B_Global_Lat   = erbx_get_bracket_info(Tstats, "Global", cfg.measureLatency);
B_Global_Mean  = erbx_get_bracket_info(Tstats, "Global", cfg.measureMean);
B_Global_Slope = erbx_get_bracket_info(Tstats, "Global", cfg.measureSlope);

fprintf(['\n[EEG RainBox] COUNTING vs VISUAL INTERFERENCE SUMMARY ' ...
    '(%s-corrected within each Effect x Measure 3-test family)\n'], ...
    upper(char(cfg.multCompMethod)));

for r = 1:height(Tstats)
    fprintf('  %s | %s | %s: n=%d, raw p=%.4f, corrected p=%.4f, %s\n', ...
        Tstats.Effect(r), Tstats.Measure(r), Tstats.Condition(r), ...
        Tstats.Npairs(r), Tstats.RawPValue(r), Tstats.CorrectedPValue(r), ...
        Tstats.Significance(r));
end

%% -------------------- figures: MEAN --------------------
figLocal = erbx_make_effect_figure( ...
    "Local", ...
    ["Local Standard","Local Deviant"], ...
    {LS_Lat_E1, LD_Lat_E1}, {LS_Lat_E2, LD_Lat_E2}, ...
    {LS_Mean_E1, LD_Mean_E1}, {LS_Mean_E2, LD_Mean_E2}, ...
    "Mean amplitude", "ERP mean (\muV)", ...
    B_Local_Lat, B_Local_Mean, cfg);

outPngLocal = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_Local_E1vsE2.png', char(cfg.chanName)));
exportgraphics(figLocal, outPngLocal, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfLocal = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_Local_E1vsE2.pdf', char(cfg.chanName)));
    exportgraphics(figLocal, outPdfLocal, 'BackgroundColor', 'white');
end
close(figLocal);

figGlobal = erbx_make_effect_figure( ...
    "Global", ...
    ["Global Standard","Global Deviant"], ...
    {GS_Lat_E1, GD_Lat_E1}, {GS_Lat_E2, GD_Lat_E2}, ...
    {GS_Mean_E1, GD_Mean_E1}, {GS_Mean_E2, GD_Mean_E2}, ...
    "Mean amplitude", "ERP mean (\muV)", ...
    B_Global_Lat, B_Global_Mean, cfg);

outPngGlobal = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_Global_E1vsE2.png', char(cfg.chanName)));
exportgraphics(figGlobal, outPngGlobal, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfGlobal = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_Global_E1vsE2.pdf', char(cfg.chanName)));
    exportgraphics(figGlobal, outPdfGlobal, 'BackgroundColor', 'white');
end
close(figGlobal);

%% -------------------- figures: SLOPE --------------------
if strcmpi(cfg.slopeUnit, "ms")
    slopeYLab = "ERP slope (\muV/ms)";
else
    slopeYLab = "ERP slope (\muV/s)";
end

figLocalSlope = erbx_make_effect_figure( ...
    "Local", ...
    ["Local Standard","Local Deviant"], ...
    {LS_Lat_E1, LD_Lat_E1}, {LS_Lat_E2, LD_Lat_E2}, ...
    {LS_Slope_E1, LD_Slope_E1}, {LS_Slope_E2, LD_Slope_E2}, ...
    "Slope", slopeYLab, ...
    B_Local_Lat, B_Local_Slope, cfg);

outPngLocalSlope = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_LocalSlope_E1vsE2.png', char(cfg.chanName)));
exportgraphics(figLocalSlope, outPngLocalSlope, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfLocalSlope = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_LocalSlope_E1vsE2.pdf', char(cfg.chanName)));
    exportgraphics(figLocalSlope, outPdfLocalSlope, 'BackgroundColor', 'white');
end
close(figLocalSlope);

figGlobalSlope = erbx_make_effect_figure( ...
    "Global", ...
    ["Global Standard","Global Deviant"], ...
    {GS_Lat_E1, GD_Lat_E1}, {GS_Lat_E2, GD_Lat_E2}, ...
    {GS_Slope_E1, GD_Slope_E1}, {GS_Slope_E2, GD_Slope_E2}, ...
    "Slope", slopeYLab, ...
    B_Global_Lat, B_Global_Slope, cfg);

outPngGlobalSlope = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_GlobalSlope_E1vsE2.png', char(cfg.chanName)));
exportgraphics(figGlobalSlope, outPngGlobalSlope, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    outPdfGlobalSlope = fullfile(grandFold, sprintf('EEG_RAINBOX_COMPARE_%s_GlobalSlope_E1vsE2.pdf', char(cfg.chanName)));
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

fprintf('\n[EEG RainBox] Saved:\n');
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

function fallbackWin = erbx_get_fallback_window(times)
times = double(times(:)');
if isempty(times)
    fallbackWin = [0 700];
else
    fallbackWin = [times(1) times(end)];
end
end

function W = erbx_mask_span_window(relTimes, sigMask, fallbackWin)
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

function subs = erbx_list_subjects_with_file(outroot, expKey, excludeSubs)
d = dir(fullfile(outroot, 'S*'));
d = d([d.isdir]);

names = string({d.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+$', 'once'));
subs = names(keep);

if ~isempty(excludeSubs)
    subs = subs(~ismember(lower(subs), lower(excludeSubs)));
end

hasFile = false(size(subs));
for i = 1:numel(subs)
    f = fullfile(outroot, char(subs(i)), sprintf('auditory_all_LG_%s_%s.set', char(subs(i)), expKey));
    hasFile(i) = exist(f,'file') == 2;
end

subs = subs(hasFile);
subs = erbx_sort_subject_ids(subs);
end

function subs = erbx_sort_subject_ids(subs)
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

function C = erbx_load_subject_curves_eeg(setFile, chanName, targets)
C = [];

EEG = pop_loadset(setFile);
EEG = eeg_checkset(EEG);

if ~(isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type'))
    return;
end

ch = find(strcmpi({EEG.chanlocs.labels}, char(chanName)), 1);
if isempty(ch)
    warning('Channel %s not found in %s. Skipping.', char(chanName), setFile);
    return;
end

M = false(4, EEG.trials);
for ep = 1:EEG.trials
    lablist = erbx_to_strlist(EEG.epoch(ep).eventlg_type);
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

FLS = squeeze(mean(double(EEG.data(ch,:,idxFLS)),3)); FLS = FLS(:)';
RLS = squeeze(mean(double(EEG.data(ch,:,idxRLS)),3)); RLS = RLS(:)';
FLD = squeeze(mean(double(EEG.data(ch,:,idxFLD)),3)); FLD = FLD(:)';
RLD = squeeze(mean(double(EEG.data(ch,:,idxRLD)),3)); RLD = RLD(:)';

C = struct();
C.relTimes = double(EEG.times(:))';

C.FLS = FLS;
C.RLS = RLS;
C.FLD = FLD;
C.RLD = RLD;

C.LS = mean([FLS; RLS], 1, 'omitnan');
C.LD = mean([FLD; RLD], 1, 'omitnan');

C.GS = mean([FLS; FLD], 1, 'omitnan');
C.GD = mean([RLS; RLD], 1, 'omitnan');
end

function S = erbx_to_strlist(v)
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

function val = erbx_mean_in_window(relTimes, curve, startMS, endMS)
idx = relTimes >= startMS & relTimes <= endMS;
if ~any(idx) || all(isnan(curve(idx)))
    val = NaN;
else
    val = mean(curve(idx), 'omitnan');
end
end

function lat = erbx_curve_peak_latency(relTimes, curve, startMS, endMS, peakPolarity)
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

switch lower(char(peakPolarity))
    case 'negative'
        [~,ix] = min(y);
    case 'absmax'
        [~,ix] = max(abs(y));
    otherwise
        [~,ix] = max(y);   % positive peak at Pz
end

lat = t(ix);
end

function slopeVal = erbx_curve_slope(relTimes, curve, startMS, endMS, slopeMode, slopeUnit, peakPolarity)
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

    otherwise  % start_to_peak
        switch lower(char(peakPolarity))
            case 'negative'
                [~,ixPeak] = min(y);
            case 'absmax'
                [~,ixPeak] = max(abs(y));
            otherwise
                [~,ixPeak] = max(y);
        end
        tFit = t(1:ixPeak);
        yFit = y(1:ixPeak);
end

if numel(tFit) < 2 || numel(unique(tFit)) < 2
    slopeVal = NaN;
    return;
end

switch lower(char(slopeUnit))
    case 'ms'
        x = tFit - tFit(1);
    otherwise
        x = (tFit - tFit(1)) / 1000;
end

p = polyfit(x, yFit, 1);
slopeVal = p(1);
end

function Twin = erbx_build_window_table(Wexp, cfg)
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

function Tstats = erbx_build_stats_table(cfg, ...
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
row=row+1; [p,~,n] = erbx_paired_pvalue(LS_Lat_E1, LS_Lat_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Peak latency"; Condition(row)="Local Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(LD_Lat_E1, LD_Lat_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Peak latency"; Condition(row)="Local Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(LD_Lat_E1 - LS_Lat_E1, LD_Lat_E2 - LS_Lat_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Peak latency"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- LOCAL MEAN ----------
row=row+1; [p,~,n] = erbx_paired_pvalue(LS_Mean_E1, LS_Mean_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Mean amplitude"; Condition(row)="Local Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(LD_Mean_E1, LD_Mean_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Mean amplitude"; Condition(row)="Local Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(LD_Mean_E1 - LS_Mean_E1, LD_Mean_E2 - LS_Mean_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Mean amplitude"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- LOCAL SLOPE ----------
row=row+1; [p,~,n] = erbx_paired_pvalue(LS_Slope_E1, LS_Slope_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Slope amplitude"; Condition(row)="Local Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(LD_Slope_E1, LD_Slope_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Slope amplitude"; Condition(row)="Local Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(LD_Slope_E1 - LS_Slope_E1, LD_Slope_E2 - LS_Slope_E2, cfg.statsMethod);
Effect(row)="Local"; Measure(row)="Slope amplitude"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- GLOBAL LATENCY ----------
row=row+1; [p,~,n] = erbx_paired_pvalue(GS_Lat_E1, GS_Lat_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Peak latency"; Condition(row)="Global Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(GD_Lat_E1, GD_Lat_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Peak latency"; Condition(row)="Global Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(GD_Lat_E1 - GS_Lat_E1, GD_Lat_E2 - GS_Lat_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Peak latency"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- GLOBAL MEAN ----------
row=row+1; [p,~,n] = erbx_paired_pvalue(GS_Mean_E1, GS_Mean_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Mean amplitude"; Condition(row)="Global Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(GD_Mean_E1, GD_Mean_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Mean amplitude"; Condition(row)="Global Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(GD_Mean_E1 - GS_Mean_E1, GD_Mean_E2 - GS_Mean_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Mean amplitude"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- GLOBAL SLOPE ----------
row=row+1; [p,~,n] = erbx_paired_pvalue(GS_Slope_E1, GS_Slope_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Slope amplitude"; Condition(row)="Global Standard"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(GD_Slope_E1, GD_Slope_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Slope amplitude"; Condition(row)="Global Deviant"; Npairs(row)=n; RawPValue(row)=p;

row=row+1; [p,~,n] = erbx_paired_pvalue(GD_Slope_E1 - GS_Slope_E1, GD_Slope_E2 - GS_Slope_E2, cfg.statsMethod);
Effect(row)="Global"; Measure(row)="Slope amplitude"; Condition(row)="Deviant - Standard contrast"; Npairs(row)=n; RawPValue(row)=p;

% ---------- family labels ----------
for k = 1:numel(Effect)
    Family(k) = Effect(k) + " | " + Measure(k);
end

% ---------- correction within each 3-test family ----------
uFam = unique(Family, 'stable');
for iFam = 1:numel(uFam)
    idxFam = (Family == uFam(iFam));
    CorrectedPValue(idxFam) = erbx_correct_pvalues(RawPValue(idxFam), cfg.multCompMethod);
end

for k = 1:numel(CorrectedPValue)
    Significance(k) = erbx_p_to_stars(CorrectedPValue(k));
end

Tstats = table( ...
    Effect, Measure, Family, Condition, ExperimentA, ExperimentB, ...
    Npairs, RawPValue, CorrectedPValue, Significance);
end

function fig = erbx_make_effect_figure(effectName, condLabels, ...
    latE1Cell, latE2Cell, rightE1Cell, rightE2Cell, ...
    rightPanelTitle, rightYLabel, BLat, BRight, cfg)

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
erbx_style_axis(ax1, cfg.bgLatency);

title(ax1, sprintf('%s effects | ERP peak latency', effectName), 'FontSize', 12, 'FontWeight','bold');
ylabel(ax1, 'ERP latency (ms)', 'FontSize', 12);

latSupport = [];
latDataAll = [];

for c = 1:2
    y1 = latE1Cell{c};
    y2 = latE2Cell{c};

    [v1min,v1max] = erbx_draw_half_violin(ax1, y1, xA(c), 'left',  cfg.colorE1, cfg.rainWidth, cfg.violinAlpha);
    [v2min,v2max] = erbx_draw_half_violin(ax1, y2, xB(c), 'right', cfg.colorE2, cfg.rainWidth, cfg.violinAlpha);

    latSupport = [latSupport; v1min; v1max; v2min; v2max]; %#ok<AGROW>
    latDataAll = [latDataAll; y1(:); y2(:)]; %#ok<AGROW>

    erbx_draw_paired_lines(ax1, y1, y2, xA(c), xB(c), cfg.grayLine);
    erbx_draw_custom_box(ax1, y1, xA(c), cfg.colorE1, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    erbx_draw_custom_box(ax1, y2, xB(c), cfg.colorE2, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    erbx_draw_points(ax1, y1, xA(c), cfg.colorE1, cfg.jitter, cfg.pointSize);
    erbx_draw_points(ax1, y2, xB(c), cfg.colorE2, cfg.jitter, cfg.pointSize);
end

set(ax1, 'XLim', [0.5 2.5], 'XTick', groupCenters, 'XTickLabel', cellstr(condLabels));
erbx_set_nice_ylim(ax1, latDataAll, latSupport);

erbx_add_pair_bracket(ax1, xA(1), xB(1), [latE1Cell{1}; latE2Cell{1}], 1, BLat.cond1);
erbx_add_pair_bracket(ax1, xA(2), xB(2), [latE1Cell{2}; latE2Cell{2}], 2, BLat.cond2);
erbx_add_global_bracket(ax1, mean([xA(1) xB(1)]), mean([xA(2) xB(2)]), latDataAll, BLat.contrast);

%% ---------- RIGHT: GENERIC PARAMETER ----------
ax2 = nexttile(tl); hold(ax2,'on');
erbx_style_axis(ax2, cfg.bgERP);

title(ax2, sprintf('%s effects | ERP %s', effectName, lower(rightPanelTitle)), 'FontSize', 12, 'FontWeight','bold');
ylabel(ax2, rightYLabel, 'FontSize', 12);

rightSupport = [];
rightDataAll = [];

for c = 1:2
    y1 = rightE1Cell{c};
    y2 = rightE2Cell{c};

    [v1min,v1max] = erbx_draw_half_violin(ax2, y1, xA(c), 'left',  cfg.colorE1, cfg.rainWidth, cfg.violinAlpha);
    [v2min,v2max] = erbx_draw_half_violin(ax2, y2, xB(c), 'right', cfg.colorE2, cfg.rainWidth, cfg.violinAlpha);

    rightSupport = [rightSupport; v1min; v1max; v2min; v2max]; %#ok<AGROW>
    rightDataAll = [rightDataAll; y1(:); y2(:)]; %#ok<AGROW>

    erbx_draw_paired_lines(ax2, y1, y2, xA(c), xB(c), cfg.grayLine);
    erbx_draw_custom_box(ax2, y1, xA(c), cfg.colorE1, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    erbx_draw_custom_box(ax2, y2, xB(c), cfg.colorE2, cfg.boxWidth, cfg.boxLineWidth, cfg.meanLineWidth);
    erbx_draw_points(ax2, y1, xA(c), cfg.colorE1, cfg.jitter, cfg.pointSize);
    erbx_draw_points(ax2, y2, xB(c), cfg.colorE2, cfg.jitter, cfg.pointSize);
end

set(ax2, 'XLim', [0.5 2.5], 'XTick', groupCenters, 'XTickLabel', cellstr(condLabels));
erbx_set_nice_ylim(ax2, rightDataAll, rightSupport);

erbx_add_pair_bracket(ax2, xA(1), xB(1), [rightE1Cell{1}; rightE2Cell{1}], 1, BRight.cond1);
erbx_add_pair_bracket(ax2, xA(2), xB(2), [rightE1Cell{2}; rightE2Cell{2}], 2, BRight.cond2);
erbx_add_global_bracket(ax2, mean([xA(1) xB(1)]), mean([xA(2) xB(2)]), rightDataAll, BRight.contrast);

%% ---------- Legend LEFT ----------
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

%% ---------- Legend RIGHT ----------
h1R = plot(ax2, nan, nan, '-', 'Color', cfg.colorE1, 'LineWidth', 3, 'DisplayName', char(cfg.labelE1));
h2R = plot(ax2, nan, nan, '-', 'Color', cfg.colorE2, 'LineWidth', 3, 'DisplayName', char(cfg.labelE2));

lgd2 = legend(ax2, [h1R h2R], 'Location', 'none', 'Box', 'off');
lgd2.FontSize = 10;
lgd2.Units = 'normalized';

ax2Pos = ax2.Position;
legendX = ax2Pos(1) + ax2Pos(3)/2 - legendW/2;
legendY = ax2Pos(2) + 0.03;
lgd2.Position = [legendX legendY legendW legendH];

sgtitle(fig, sprintf('Event Related Potentials at %s | %s effects', upper(char(cfg.chanName)), effectName), ...
    'FontSize', 12, 'FontWeight', 'bold');

end

function erbx_style_axis(ax, bgColor)
set(ax, ...
    'Color', bgColor, ...
    'FontSize', 11, ...
    'LineWidth', 1, ...
    'Box', 'off', ...
    'Layer', 'top', ...
    'TickDir', 'out');
grid(ax, 'off');
end

function [vmin, vmax] = erbx_draw_half_violin(ax, y, x0, sideStr, rgb, maxWidth, faceAlpha)
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

function erbx_draw_paired_lines(ax, y1, y2, x1, x2, lineColor)
m = isfinite(y1) & isfinite(y2);
idx = find(m);
for i = 1:numel(idx)
    ii = idx(i);
    plot(ax, [x1 x2], [y1(ii) y2(ii)], '-', ...
        'Color', lineColor, 'LineWidth', 1.5, 'HandleVisibility','off');
end
end

function erbx_draw_points(ax, y, x0, rgb, jitterAmt, markerSize)
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

function erbx_draw_custom_box(ax, y, x0, rgb, boxW, boxLW, meanLW)
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

function erbx_set_nice_ylim(ax, yData, ySupport)
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

function yTop = erbx_add_pair_bracket(ax, x1, x2, yVals, groupIdx, txt)
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

function erbx_add_global_bracket(ax, x1, x2, ~, txt)
yl = ylim(ax);
yr = yl(2) - yl(1);

y = yl(2) - 0.10 * yr; % second number, the smaller the closer to top
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

function [p, stars, n] = erbx_paired_pvalue(a, b, methodStr)
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

stars = erbx_p_to_stars(p);
end

function txt = erbx_p_to_stars(p)
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

function B = erbx_get_bracket_info(Tstats, effectName, measureName)
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

function pAdj = erbx_correct_pvalues(pRaw, methodStr)
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

    otherwise % holm
        pAdj(ok) = erbx_holm_adjust(p);
end
end

function pAdj = erbx_holm_adjust(p)
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