function OUT = pupil_export_rainbox_compare_E1E2_ANOVA(OUT_RAINBOX, grandOutFolder)
% pupil_export_rainbox_compare_E1E2_ANOVA
%
% Separate 2 x 2 repeated-measures ANOVA analysis for the pupil rainbox
% summaries. This function does not modify the original Wilcoxon function.
%
% INPUTS
%   OUT_RAINBOX      Output returned by pupil_export_rainbox_compare_E1E2
%   grandOutFolder   Folder in which figures and result tables are saved
%
% DESIGN (run separately for Local/Global and for each measure)
%   Within-subject factor 1: Task      (Counting, Visual interference)
%   Within-subject factor 2: Condition (Standard, Deviant)
%
% FIGURE BRACKETS
%   Short brackets: paired Task simple effects for Standard and Deviant.
%   Long bracket:   Task x Condition interaction from the 2 x 2 ANOVA.
%   The three displayed p-values are Holm-adjusted together within each
%   Effect x Measure family, matching the original three-bracket logic.
%
% OUTPUTS
%   Four figures with the same visual style as the original rainbox plots,
%   but with filenames containing "ANOVA".
%   CSV/MAT files containing omnibus ANOVA and follow-up results.
%
% Requires Statistics and Machine Learning Toolbox (fitrm, ranova, ttest).
%
% Miguel Contreras-Altamirano, 2026

%% -------------------- input checks --------------------
if nargin < 2 || isempty(grandOutFolder)
    error('Both OUT_RAINBOX and grandOutFolder are required.');
end
if ~isstruct(OUT_RAINBOX) || ~isfield(OUT_RAINBOX,'subjectTable')
    error('OUT_RAINBOX must contain the field subjectTable.');
end
if ~exist(grandOutFolder,'dir')
    mkdir(grandOutFolder);
end

Tsub = OUT_RAINBOX.subjectTable;

%% -------------------- style/config (same as original figure) --------------------
cfg = struct();
cfg.labelE1 = "Counting";
cfg.labelE2 = "Visual interference";

cfg.colorE1 = [0 0.45 0.85];
cfg.colorE2 = [0.90 0.45 0.10];
cfg.grayLine = [0.62 0.62 0.62];

cfg.bgLatency = [1 1 1];
cfg.bgPupil   = [1 1 1];

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

cfg.legendLocation = 'southwest';
cfg.alpha = 0.05;
cfg.followupCorrection = "holm";

% Keep labels consistent with the original function.
cfg.slopeLabel = "Pupil slope (z-score/ms)";
if isfield(OUT_RAINBOX,'cfg') && isfield(OUT_RAINBOX.cfg,'slopeLabel')
    cfg.slopeLabel = string(OUT_RAINBOX.cfg.slopeLabel);
end

fprintf('\n============================================================\n');
fprintf('PUPIL 2 x 2 REPEATED-MEASURES ANOVA (RAINBOX SUMMARIES)\n');
fprintf('Task: Counting vs Visual interference\n');
fprintf('Condition: Standard vs Deviant\n');
fprintf('Figure brackets: Holm-adjusted across Standard, Deviant, and interaction within each family\n');
fprintf('============================================================\n');

%% -------------------- collect result tables --------------------
Tanova = table();
Tfollow = table();
Tdisplay = table();
RANOVA = struct();

% LOCAL | PEAK LATENCY
[A,F,D,B,R] = anv_run_family( ...
    "Local", "Peak latency", ...
    Tsub.LocalStandard_Latency_E1_ms, Tsub.LocalDeviant_Latency_E1_ms, ...
    Tsub.LocalStandard_Latency_E2_ms, Tsub.LocalDeviant_Latency_E2_ms, ...
    cfg);
Tanova = A; Tfollow = F; Tdisplay = D;
B_Local_Lat = B; RANOVA.Local_PeakLatency = R;

% LOCAL | MEAN PUPIL SIZE
[A,F,D,B,R] = anv_run_family( ...
    "Local", "Mean pupil size", ...
    Tsub.LocalStandard_Mean_E1, Tsub.LocalDeviant_Mean_E1, ...
    Tsub.LocalStandard_Mean_E2, Tsub.LocalDeviant_Mean_E2, ...
    cfg);
Tanova = [Tanova; A]; Tfollow = [Tfollow; F]; Tdisplay = [Tdisplay; D];
B_Local_Mean = B; RANOVA.Local_MeanPupilSize = R;

% LOCAL | SLOPE
[A,F,D,B,R] = anv_run_family( ...
    "Local", "Slope pupil size", ...
    Tsub.LocalStandard_Slope_E1, Tsub.LocalDeviant_Slope_E1, ...
    Tsub.LocalStandard_Slope_E2, Tsub.LocalDeviant_Slope_E2, ...
    cfg);
Tanova = [Tanova; A]; Tfollow = [Tfollow; F]; Tdisplay = [Tdisplay; D];
B_Local_Slope = B; RANOVA.Local_SlopePupilSize = R;

% GLOBAL | PEAK LATENCY
[A,F,D,B,R] = anv_run_family( ...
    "Global", "Peak latency", ...
    Tsub.GlobalStandard_Latency_E1_ms, Tsub.GlobalDeviant_Latency_E1_ms, ...
    Tsub.GlobalStandard_Latency_E2_ms, Tsub.GlobalDeviant_Latency_E2_ms, ...
    cfg);
Tanova = [Tanova; A]; Tfollow = [Tfollow; F]; Tdisplay = [Tdisplay; D];
B_Global_Lat = B; RANOVA.Global_PeakLatency = R;

% GLOBAL | MEAN PUPIL SIZE
[A,F,D,B,R] = anv_run_family( ...
    "Global", "Mean pupil size", ...
    Tsub.GlobalStandard_Mean_E1, Tsub.GlobalDeviant_Mean_E1, ...
    Tsub.GlobalStandard_Mean_E2, Tsub.GlobalDeviant_Mean_E2, ...
    cfg);
Tanova = [Tanova; A]; Tfollow = [Tfollow; F]; Tdisplay = [Tdisplay; D];
B_Global_Mean = B; RANOVA.Global_MeanPupilSize = R;

% GLOBAL | SLOPE
[A,F,D,B,R] = anv_run_family( ...
    "Global", "Slope pupil size", ...
    Tsub.GlobalStandard_Slope_E1, Tsub.GlobalDeviant_Slope_E1, ...
    Tsub.GlobalStandard_Slope_E2, Tsub.GlobalDeviant_Slope_E2, ...
    cfg);
Tanova = [Tanova; A]; Tfollow = [Tfollow; F]; Tdisplay = [Tdisplay; D];
B_Global_Slope = B; RANOVA.Global_SlopePupilSize = R;

%% -------------------- save result tables --------------------
outAnovaCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_omnibus_results.csv');
outFollowCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_followup_results.csv');
outDisplayCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_display_pvalues.csv');
outSubjectCsv = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_subject_values.csv');
outMat = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_results.mat');

writetable(Tanova, outAnovaCsv);
writetable(Tfollow, outFollowCsv);
writetable(Tdisplay, outDisplayCsv);
writetable(Tsub, outSubjectCsv);

save(outMat, 'Tanova', 'Tfollow', 'Tdisplay', 'Tsub', 'RANOVA', 'cfg', '-v7.3');

%% -------------------- print results --------------------
fprintf('\n=== OMNIBUS 2 x 2 REPEATED-MEASURES ANOVA RESULTS ===\n');
for r = 1:height(Tanova)
    fprintf('%s | %s | %s: F(%g,%g)=%.4f, p=%.6f, partial eta^2=%.4f, n=%d\n', ...
        Tanova.Effect(r), Tanova.Measure(r), Tanova.ANOVAEffect(r), ...
        Tanova.DF1(r), Tanova.DF2(r), Tanova.FValue(r), ...
        Tanova.PValue(r), Tanova.PartialEtaSquared(r), Tanova.N(r));
end

fprintf('\n=== HOLM-ADJUSTED TASK SIMPLE EFFECTS ===\n');
for r = 1:height(Tfollow)
    fprintf('%s | %s | %s: E1=%.4f +/- %.4f, E2=%.4f +/- %.4f, ', ...
        Tfollow.Effect(r), Tfollow.Measure(r), Tfollow.Condition(r), ...
        Tfollow.MeanE1(r), Tfollow.SEME1(r), ...
        Tfollow.MeanE2(r), Tfollow.SEME2(r));
    fprintf('t(%g)=%.4f, raw p=%.6f, Holm p=%.6f, dz=%.4f, %s\n', ...
        Tfollow.DF(r), Tfollow.TValue(r), Tfollow.RawPValue(r), ...
        Tfollow.HolmPValue(r), Tfollow.Dz(r), Tfollow.Significance(r));
end

fprintf('\n=== P-VALUES USED FOR FIGURE BRACKETS ===\n');
disp(Tdisplay);

%% -------------------- figures: MEAN --------------------
figLocal = rbx_make_effect_figure( ...
    "Local", ...
    ["Local Standard","Local Deviant"], ...
    {Tsub.LocalStandard_Latency_E1_ms, Tsub.LocalDeviant_Latency_E1_ms}, ...
    {Tsub.LocalStandard_Latency_E2_ms, Tsub.LocalDeviant_Latency_E2_ms}, ...
    {Tsub.LocalStandard_Mean_E1, Tsub.LocalDeviant_Mean_E1}, ...
    {Tsub.LocalStandard_Mean_E2, Tsub.LocalDeviant_Mean_E2}, ...
    "mean", "Pupil mean (z-score)", ...
    B_Local_Lat, B_Local_Mean, cfg);

outPngLocal = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_Local_E1vsE2.png');
exportgraphics(figLocal, outPngLocal, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    exportgraphics(figLocal, replace(outPngLocal,'.png','.pdf'), 'BackgroundColor', 'white');
end
close(figLocal);

figGlobal = rbx_make_effect_figure( ...
    "Global", ...
    ["Global Standard","Global Deviant"], ...
    {Tsub.GlobalStandard_Latency_E1_ms, Tsub.GlobalDeviant_Latency_E1_ms}, ...
    {Tsub.GlobalStandard_Latency_E2_ms, Tsub.GlobalDeviant_Latency_E2_ms}, ...
    {Tsub.GlobalStandard_Mean_E1, Tsub.GlobalDeviant_Mean_E1}, ...
    {Tsub.GlobalStandard_Mean_E2, Tsub.GlobalDeviant_Mean_E2}, ...
    "mean", "Pupil mean (z-score)", ...
    B_Global_Lat, B_Global_Mean, cfg);

outPngGlobal = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_Global_E1vsE2.png');
exportgraphics(figGlobal, outPngGlobal, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    exportgraphics(figGlobal, replace(outPngGlobal,'.png','.pdf'), 'BackgroundColor', 'white');
end
close(figGlobal);

%% -------------------- figures: SLOPE --------------------
figLocalSlope = rbx_make_effect_figure( ...
    "Local", ...
    ["Local Standard","Local Deviant"], ...
    {Tsub.LocalStandard_Latency_E1_ms, Tsub.LocalDeviant_Latency_E1_ms}, ...
    {Tsub.LocalStandard_Latency_E2_ms, Tsub.LocalDeviant_Latency_E2_ms}, ...
    {Tsub.LocalStandard_Slope_E1, Tsub.LocalDeviant_Slope_E1}, ...
    {Tsub.LocalStandard_Slope_E2, Tsub.LocalDeviant_Slope_E2}, ...
    "Slope", char(cfg.slopeLabel), ...
    B_Local_Lat, B_Local_Slope, cfg);

outPngLocalSlope = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_LocalSlope_E1vsE2.png');
exportgraphics(figLocalSlope, outPngLocalSlope, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    exportgraphics(figLocalSlope, replace(outPngLocalSlope,'.png','.pdf'), 'BackgroundColor', 'white');
end
close(figLocalSlope);

figGlobalSlope = rbx_make_effect_figure( ...
    "Global", ...
    ["Global Standard","Global Deviant"], ...
    {Tsub.GlobalStandard_Latency_E1_ms, Tsub.GlobalDeviant_Latency_E1_ms}, ...
    {Tsub.GlobalStandard_Latency_E2_ms, Tsub.GlobalDeviant_Latency_E2_ms}, ...
    {Tsub.GlobalStandard_Slope_E1, Tsub.GlobalDeviant_Slope_E1}, ...
    {Tsub.GlobalStandard_Slope_E2, Tsub.GlobalDeviant_Slope_E2}, ...
    "Slope", char(cfg.slopeLabel), ...
    B_Global_Lat, B_Global_Slope, cfg);

outPngGlobalSlope = fullfile(grandOutFolder, 'PUPIL_RAINBOX_COMPARE_ANOVA_GlobalSlope_E1vsE2.png');
exportgraphics(figGlobalSlope, outPngGlobalSlope, 'Resolution', cfg.dpi, 'BackgroundColor', 'white');
if cfg.savePDF
    exportgraphics(figGlobalSlope, replace(outPngGlobalSlope,'.png','.pdf'), 'BackgroundColor', 'white');
end
close(figGlobalSlope);

%% -------------------- output structure --------------------
OUT = struct();
OUT.subjectTable = Tsub;
OUT.anovaTable = Tanova;
OUT.followupTable = Tfollow;
OUT.displayPTable = Tdisplay;
OUT.ranovaTables = RANOVA;
OUT.cfg = cfg;

if isfield(OUT_RAINBOX,'windowTable')
    OUT.windowTable = OUT_RAINBOX.windowTable;
end
if isfield(OUT_RAINBOX,'windowsExperiment')
    OUT.windowsExperiment = OUT_RAINBOX.windowsExperiment;
end

OUT.files.localPNG       = outPngLocal;
OUT.files.globalPNG      = outPngGlobal;
OUT.files.localSlopePNG  = outPngLocalSlope;
OUT.files.globalSlopePNG = outPngGlobalSlope;
OUT.files.anovaCSV       = outAnovaCsv;
OUT.files.followupCSV    = outFollowCsv;
OUT.files.displayCSV     = outDisplayCsv;
OUT.files.subjectCSV     = outSubjectCsv;
OUT.files.resultsMAT     = outMat;

fprintf('\n[ANOVA RainBox] Saved:\n');
fprintf('  %s\n', outPngLocal);
fprintf('  %s\n', outPngGlobal);
fprintf('  %s\n', outPngLocalSlope);
fprintf('  %s\n', outPngGlobalSlope);
fprintf('  %s\n', outAnovaCsv);
fprintf('  %s\n', outFollowCsv);
fprintf('  %s\n', outDisplayCsv);
fprintf('  %s\n', outMat);
end

%% ========================================================================
% ANOVA HELPERS
% ========================================================================

function [Tanova, Tfollow, Tdisplay, B, R] = anv_run_family(effectName, measureName, stdE1, devE1, stdE2, devE2, cfg)

Y = [double(stdE1(:)), double(devE1(:)), double(stdE2(:)), double(devE2(:))];
keep = all(isfinite(Y),2);
Y = Y(keep,:);
n = size(Y,1);

if n < 3
    error('Insufficient complete cases for %s | %s (n=%d).', effectName, measureName, n);
end

% Four repeated measures in the order:
% E1 Standard, E1 Deviant, E2 Standard, E2 Deviant
Twide = array2table(Y, 'VariableNames', {'Y1','Y2','Y3','Y4'});
Within = table( ...
    categorical(["Counting"; "Counting"; "VisualInterference"; "VisualInterference"]), ...
    categorical(["Standard"; "Deviant"; "Standard"; "Deviant"]), ...
    'VariableNames', {'Task','Condition'});

% Actual 2 x 2 repeated-measures ANOVA.
R = table();
usedFitrm = true;
try
    rm = fitrm(Twide, 'Y1-Y4 ~ 1', 'WithinDesign', Within);
    R = ranova(rm, 'WithinModel', 'Task*Condition');
catch
    usedFitrm = false;
    warning('fitrm/ranova failed for %s | %s. Using exact 1-df contrast form.', ...
        effectName, measureName);
end

% Direct contrast formulation. For a 2 x 2 within-subject design, each
% effect has 1 numerator df and F = t^2. This is also used as a robust
% fallback and as an internal check on fitrm/ranova.
Ctask = mean(Y(:,1:2),2) - mean(Y(:,3:4),2);
Ccond = mean(Y(:,[2 4]),2) - mean(Y(:,[1 3]),2);
Cint  = (Y(:,2)-Y(:,1)) - (Y(:,4)-Y(:,3));

[Ftask, pTask, dfTask] = anv_contrast_F(Ctask);
[Fcond, pCond, dfCond] = anv_contrast_F(Ccond);
[Fint,  pInt,  dfInt]  = anv_contrast_F(Cint);

if usedFitrm
    [FtaskR,pTaskR] = anv_extract_ranova_effect(R,"Task");
    [FcondR,pCondR] = anv_extract_ranova_effect(R,"Condition");
    [FintR,pIntR]   = anv_extract_ranova_effect(R,"Interaction");

    if isfinite(FtaskR), Ftask = FtaskR; end
    if isfinite(pTaskR), pTask = pTaskR; end
    if isfinite(FcondR), Fcond = FcondR; end
    if isfinite(pCondR), pCond = pCondR; end
    if isfinite(FintR),  Fint  = FintR;  end
    if isfinite(pIntR),  pInt  = pIntR;  end
end

df1 = 1;
df2 = n - 1;
etaTask = anv_partial_eta_squared(Ftask,df1,df2);
etaCond = anv_partial_eta_squared(Fcond,df1,df2);
etaInt  = anv_partial_eta_squared(Fint, df1,df2);

Tanova = table( ...
    repmat(string(effectName),3,1), ...
    repmat(string(measureName),3,1), ...
    ["Task"; "Condition"; "Task x Condition"], ...
    repmat(n,3,1), ...
    repmat(df1,3,1), ...
    repmat(df2,3,1), ...
    [Ftask; Fcond; Fint], ...
    [pTask; pCond; pInt], ...
    [etaTask; etaCond; etaInt], ...
    'VariableNames', {'Effect','Measure','ANOVAEffect','N','DF1','DF2', ...
    'FValue','PValue','PartialEtaSquared'});

% Follow-up paired Task comparisons at each Condition level.
[s1, p1, t1, df1s, dz1, m11, se11, m12, se12] = anv_simple_effect(Y(:,1),Y(:,3));
[s2, p2, t2, df2s, dz2, m21, se21, m22, se22] = anv_simple_effect(Y(:,2),Y(:,4));

% Match the original three-bracket family: Standard, Deviant, and the
% Task x Condition interaction are Holm-adjusted together.
pAdjDisplay = rbx_holm_adjust([p1; p2; pInt]);

if effectName == "Local"
    condNames = ["Local Standard"; "Local Deviant"];
else
    condNames = ["Global Standard"; "Global Deviant"];
end

Tfollow = table( ...
    repmat(string(effectName),2,1), ...
    repmat(string(measureName),2,1), ...
    condNames, ...
    [s1; s2], ...
    [m11; m21], [se11; se21], [m12; m22], [se12; se22], ...
    [t1; t2], [df1s; df2s], [p1; p2], pAdjDisplay(1:2), [dz1; dz2], ...
    [rbx_p_to_stars(pAdjDisplay(1)); rbx_p_to_stars(pAdjDisplay(2))], ...
    'VariableNames', {'Effect','Measure','Condition','N', ...
    'MeanE1','SEME1','MeanE2','SEME2','TValue','DF', ...
    'RawPValue','HolmPValue','Dz','Significance'});

% The figure keeps exactly the same three brackets as the original:
% Standard simple effect, Deviant simple effect, and Task x Condition.
Tdisplay = table( ...
    repmat(string(effectName),3,1), ...
    repmat(string(measureName),3,1), ...
    [condNames; "Task x Condition interaction"], ...
    [p1; p2; pInt], pAdjDisplay, ...
    ["Paired t-test"; "Paired t-test"; "2 x 2 repeated-measures ANOVA"], ...
    [rbx_p_to_stars(pAdjDisplay(1)); rbx_p_to_stars(pAdjDisplay(2)); rbx_p_to_stars(pAdjDisplay(3))], ...
    'VariableNames', {'Effect','Measure','Bracket','RawPValue','HolmPValue','Source','Significance'});

B = struct();
B.cond1 = rbx_p_to_stars(pAdjDisplay(1));
B.cond2 = rbx_p_to_stars(pAdjDisplay(2));
B.contrast = rbx_p_to_stars(pAdjDisplay(3));
B.p1 = pAdjDisplay(1);
B.p2 = pAdjDisplay(2);
B.p3 = pAdjDisplay(3);
end

function [F,p,df] = anv_contrast_F(x)
x = x(isfinite(x));
if numel(x) < 2
    F = NaN; p = NaN; df = NaN;
    return;
end
if all(abs(x) < 1e-12)
    F = 0; p = 1; df = numel(x)-1;
    return;
end
[~,p,~,st] = ttest(x,0);
F = st.tstat.^2;
df = st.df;
end

function [F,p] = anv_extract_ranova_effect(R,effectName)
F = NaN; p = NaN;
if isempty(R) || height(R)==0 || isempty(R.Properties.RowNames)
    return;
end

rn = lower(string(R.Properties.RowNames));
isError = contains(rn,"error");
hasTask = contains(rn,"task");
hasCond = contains(rn,"condition");

switch lower(char(effectName))
    case 'task'
        idx = find(~isError & hasTask & ~hasCond,1,'first');
    case 'condition'
        idx = find(~isError & hasCond & ~hasTask,1,'first');
    otherwise
        idx = find(~isError & hasTask & hasCond,1,'first');
end

if isempty(idx)
    return;
end

if ismember('F',R.Properties.VariableNames)
    F = double(R.F(idx));
end

if ismember('pValue',R.Properties.VariableNames)
    p = double(R.pValue(idx));
elseif ismember('pValueGG',R.Properties.VariableNames)
    p = double(R.pValueGG(idx));
end
end

function eta = anv_partial_eta_squared(F,df1,df2)
if ~isfinite(F) || ~isfinite(df1) || ~isfinite(df2)
    eta = NaN;
else
    eta = (F*df1) / (F*df1 + df2);
end
end

function [n,p,tval,df,dz,m1,se1,m2,se2] = anv_simple_effect(a,b)
a = double(a(:));
b = double(b(:));
keep = isfinite(a) & isfinite(b);
a = a(keep); b = b(keep);
n = numel(a);

m1 = mean(a,'omitnan');
m2 = mean(b,'omitnan');
se1 = std(a,0,'omitnan')/sqrt(n);
se2 = std(b,0,'omitnan')/sqrt(n);

if n < 2
    p=NaN; tval=NaN; df=NaN; dz=NaN;
    return;
end

if all(abs(a-b) < 1e-12)
    p=1; tval=0; df=n-1; dz=0;
    return;
end

[~,p,~,st] = ttest(a,b);
tval = st.tstat;
df = st.df;
d = a-b;
sdD = std(d,0);
if sdD > 0
    dz = mean(d)/sdD;
else
    dz = NaN;
end
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