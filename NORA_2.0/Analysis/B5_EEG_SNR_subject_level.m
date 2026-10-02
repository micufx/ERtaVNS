clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% Calculate_SNR_single_subject_dashboards_TLG.m
% Single-subject dashboards (matches FINAL group-level logic):
%   - Baseline window (noise) from T_LG.BaseWin_E1/E2 (per subject & experiment)
%   - Pooled-condition SNR: noise = SD of ALL baseline samples pooled across selected trials
%   - Effect SNR: noise = sqrt(noiseDev^2 + noiseStd^2)
%   - Plateau method selectable: 'band' or 'slope'
%   - Plot shading (green BL + blue ERP) uses display baseline window baseWin_this (for consistent visuals)
%   - 8 figures/subject: E1(Fig1..4) + E2(Fig1..4)
%   - Writes plateau-N into T_LG (8 columns): E1/E2 × {LS,LD,GS,GD}
%
% Miguel Contreras-Altamirano, 2026

%% ===== USER SETTINGS (keep consistent with GROUP script) =====
mainpath = 'L:\Downloads\eeglab2023.0\';
outpath  = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';

% ===================== NEW: RUN CONTROL =====================
startSubIndex = 30;   % <-- Set to 24 to start from the 24th S* folder (after sorting)
% ============================================================

xlim_ms      = [-2500 700];       % plot x-range (ERP)
baseWin_this = [-2500 -1700];     % BASELINE shading window (display only)
dpi          = 200;

% ERP windows (relative to tone5=0, like your Bekinschtein plots)
MMN.sigWin = [120 220];   MMN.pol = 'neg';   MMN.chan = 'Fz';
P3B.sigWin = [300 600];   P3B.pol = 'pos';   P3B.chan = 'Pz';

% Bootstrap curve settings
nBoot = 100;
minN  = 10;
stepN = 5;                    % group is large; use 100 or 200

% Color plots
COL.baseRGB  = [0.5 1 0.5];       % baseline color (green)
COL.sigRGB   = [1.0 0.5 0.0];     % ERP color (blue)
COL.alphaVal = 0.10;              % transparency

% ---------- Plateau criterion ----------
plateauMethod  = 'slope';   % 'band' or 'slope'

% Band plateau params:
plateauTolBand = 0.10;
finalK         = 3;
plateauMinPts  = 3;

% Slope plateau params:
slopeTol       = 0.05;      % 5% of target per step
slopeMinPts    = 3;
slopeFinalK    = 3;

% N100 relative to subject X1 (from T_LG baseEnd)
N100.searchRelMS = [60 140];
N100.halfWidthMS = 20;
N100.pol         = 'neg';

% Output folder
dashSubFolderName = 'SNR_subject_level';  % you can rename if you want

% T_LG table
infoFile = fullfile(outpath,'Info_EEG_LG.mat');
assert(exist(infoFile,'file')==2, 'Missing: %s', infoFile);
Sinfo = load(infoFile);
assert(isfield(Sinfo,'T_LG'), 'Info_EEG_LG.mat does not contain T_LG');
T_LG = Sinfo.T_LG;

%% ===== EEGLAB =====
if exist('eeglab','file') ~= 2
    addpath(mainpath); % ONLY root
end
try
    eeglab nogui;
catch
    eeglab;
end
set(groot,'DefaultFigureVisible','off');

%% ===== ensure plateau columns exist in T_LG =====
T_LG = ensure_plateau_cols(T_LG);

%% ===== SUBJECT FOLDERS =====
subDirs = dir(fullfile(outpath,'S*'));
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name},{'.','..'}));
assert(~isempty(subDirs), 'No S* subject folders found in outpath.');

% ===================== NEW: Stable ordering + start index guard =====================
[~,ord] = sort({subDirs.name});   % stable: S01,S02,... (string sort)
subDirs = subDirs(ord);

if startSubIndex < 1 || startSubIndex > numel(subDirs)
    error('startSubIndex=%d is out of range. Must be between 1 and %d.', startSubIndex, numel(subDirs));
end
% ===================================================================================

%% ===== MAIN LOOP =====
for sd = startSubIndex:numel(subDirs)
    subName = subDirs(sd).name;
    subFold = fullfile(subDirs(sd).folder, subDirs(sd).name);

    fprintf('\n=== SINGLE SUBJECT | %s ===\n', subName);

    dashOut = fullfile(subFold, dashSubFolderName);
    if ~exist(dashOut,'dir'), mkdir(dashOut); end

    for expTagCell = {'E1','E2'}
        expTag = expTagCell{1};

        % baseline window from T_LG (USED FOR SNR NOISE)
        [bw, x1_ms, baseLen_ms] = get_basewin_from_TLG(T_LG, subName, expTag);
        if any(isnan([bw(:); x1_ms; baseLen_ms]))
            fprintf('  skip %s %s (no T_LG BaseWin)\n', subName, expTag);
            % clear plateau entries for this exp (optional)
            T_LG = clear_plateau_this_exp(T_LG, subName, expTag);
            continue;
        end

        % load EEG
        EEG = load_or_merge_LG_set(subFold, subName, expTag);
        if isempty(EEG) || EEG.trials <= 1
            fprintf('  skip %s %s (missing/empty EEG)\n', subName, expTag);
            T_LG = clear_plateau_this_exp(T_LG, subName, expTag);
            continue;
        end

        tms = EEG.times(:);

        % labels
        [locLab, globLab] = get_epoch_local_global_labels(EEG);
        idxLocalStd  = find(strcmpi(locLab,  'Local Standard'));
        idxLocalDev  = find(strcmpi(locLab,  'Local Deviant'));
        idxGlobalStd = find(strcmpi(globLab, 'Global Standard'));
        idxGlobalDev = find(strcmpi(globLab, 'Global Deviant'));

        fprintf('  [%s] trials: LS=%d LD=%d GS=%d GD=%d\n', expTag, ...
            numel(idxLocalStd), numel(idxLocalDev), numel(idxGlobalStd), numel(idxGlobalDev));

        % channels
        chF = pick_chan(EEG, 'Fz');
        chP = pick_chan(EEG, 'Pz');
        if isempty(chF) || isempty(chP)
            fprintf('  skip %s %s (missing Fz/Pz fallbacks)\n', subName, expTag);
            T_LG = clear_plateau_this_exp(T_LG, subName, expTag);
            continue;
        end

        % extract matrices
        Lstd = pick_trials(EEG.data, chF, idxLocalStd);
        Ldev = pick_trials(EEG.data, chF, idxLocalDev);
        Gstd = pick_trials(EEG.data, chP, idxGlobalStd);
        Gdev = pick_trials(EEG.data, chP, idxGlobalDev);

        % baseline stats (subject-specific bw)
        Lstd_stat = []; Ldev_stat = []; Gstd_stat = []; Gdev_stat = [];
        if ~isempty(Lstd), Lstd_stat = precompute_base_stats_fixedwin(Lstd, tms, bw); end
        if ~isempty(Ldev), Ldev_stat = precompute_base_stats_fixedwin(Ldev, tms, bw); end
        if ~isempty(Gstd), Gstd_stat = precompute_base_stats_fixedwin(Gstd, tms, bw); end
        if ~isempty(Gdev), Gdev_stat = precompute_base_stats_fixedwin(Gdev, tms, bw); end

        % ---- N100 window (needs Lstd present) ----
        sigWin_N100 = [NaN NaN];
        peakLat = NaN;
        if ~isempty(Lstd)
            erpN100 = mean(Lstd, 2, 'omitnan');
            searchWin = x1_ms + N100.searchRelMS;
            [sigWin_N100, peakLat] = pick_peak_window(erpN100, tms, searchWin, N100.halfWidthMS, 'neg');
        end

        fprintf('  [%s] BaseWin=[%.1f %.1f] (BL≈%.0f ms) | N100 peak %.1f ms | ERP=[%.1f %.1f]\n', ...
            expTag, bw(1), bw(2), baseLen_ms, peakLat, sigWin_N100(1), sigWin_N100(2));

        %% ===== FIG1: Local (Fz, MMN) =====
        fig1 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

        subplot(2,2,1);
        plot_buildup_from_matrix(Lstd, tms, xlim_ms, baseWin_this, baseLen_ms, MMN.sigWin, ...
            sprintf('%s_%s | Local Standard | ERP (Fz)', subName, expTag), COL);

        subplot(2,2,2);
        [NsLS, mLS] = plot_snr_curve_fast_safe(Lstd, tms, MMN.sigWin, MMN.pol, nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Lstd_stat);
        title(sprintf('%s_%s | Local Standard | SNR (MMN)', subName, expTag),'Interpreter','none');

        subplot(2,2,3);
        plot_buildup_from_matrix(Ldev, tms, xlim_ms, baseWin_this, baseLen_ms, MMN.sigWin, ...
            sprintf('%s_%s | Local Deviant | ERP (Fz)', subName, expTag), COL);

        subplot(2,2,4);
        [NsLD, mLD] = plot_snr_curve_fast_safe(Ldev, tms, MMN.sigWin, MMN.pol, nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Ldev_stat);
        title(sprintf('%s_%s | Local Deviant | SNR (MMN)', subName, expTag),'Interpreter','none');

        safe_export(fig1, fullfile(dashOut, sprintf('%s_%s_FIG1_LOCAL_MMN.png', subName, expTag)), dpi);
        close(fig1);

        %% ===== FIG2: Global (Pz, P3b) =====
        fig2 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

        subplot(2,2,1);
        plot_buildup_from_matrix(Gstd, tms, xlim_ms, baseWin_this, baseLen_ms, P3B.sigWin, ...
            sprintf('%s_%s | Global Standard | ERP (Pz)', subName, expTag), COL);

        subplot(2,2,2);
        [NsGS, mGS] = plot_snr_curve_fast_safe(Gstd, tms, P3B.sigWin, P3B.pol, nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Gstd_stat);
        title(sprintf('%s_%s | Global Standard | SNR (P3b)', subName, expTag),'Interpreter','none');

        subplot(2,2,3);
        plot_buildup_from_matrix(Gdev, tms, xlim_ms, baseWin_this, baseLen_ms, P3B.sigWin, ...
            sprintf('%s_%s | Global Deviant | ERP (Pz)', subName, expTag), COL);

        subplot(2,2,4);
        [NsGD, mGD] = plot_snr_curve_fast_safe(Gdev, tms, P3B.sigWin, P3B.pol, nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Gdev_stat);
        title(sprintf('%s_%s | Global Deviant | SNR (P3b)', subName, expTag),'Interpreter','none');

        safe_export(fig2, fullfile(dashOut, sprintf('%s_%s_FIG2_GLOBAL_P3B.png', subName, expTag)), dpi);
        close(fig2);

        %% ===== FIG3: Effects (Dev-Std) =====
        fig3 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

        subplot(2,2,1);
        plot_buildup_effect__SAFE(Lstd, Ldev, tms, xlim_ms, baseWin_this, baseLen_ms, MMN.sigWin, ...
            sprintf('%s_%s | LOCAL EFFECT (Dev-Std) | ERP (Fz)', subName, expTag), COL);

        subplot(2,2,2);
        plot_snr_curve_diff_fast_safe(Lstd, Ldev, tms, MMN.sigWin, MMN.pol, nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Lstd_stat, Ldev_stat);
        title(sprintf('%s_%s | LOCAL EFFECT (Dev-Std) | SNR (MMN)', subName, expTag),'Interpreter','none');

        subplot(2,2,3);
        plot_buildup_effect__SAFE(Gstd, Gdev, tms, xlim_ms, baseWin_this, baseLen_ms, P3B.sigWin, ...
            sprintf('%s_%s | GLOBAL EFFECT (Dev-Std) | ERP (Pz)', subName, expTag), COL);

        subplot(2,2,4);
        plot_snr_curve_diff_fast_safe(Gstd, Gdev, tms, P3B.sigWin, P3B.pol, nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Gstd_stat, Gdev_stat);
        title(sprintf('%s_%s | GLOBAL EFFECT (Dev-Std) | SNR (P3b)', subName, expTag),'Interpreter','none');

        safe_export(fig3, fullfile(dashOut, sprintf('%s_%s_FIG3_EFFECTS.png', subName, expTag)), dpi);
        close(fig3);

        %% ===== FIG4: N100 =====
        fig4 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

        subplot(2,2,1);
        plot_buildup_from_matrix(Lstd, tms, xlim_ms, baseWin_this, baseLen_ms, sigWin_N100, ...
            sprintf('%s_%s | Local Standard | ERP (Fz, N100)', subName, expTag), COL);

        subplot(2,2,2);
        plot_snr_curve_fast_safe(Lstd, tms, sigWin_N100, 'neg', nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Lstd_stat);
        title(sprintf('%s_%s | Local Standard | SNR (N100)', subName, expTag),'Interpreter','none');

        subplot(2,2,3);
        plot_buildup_from_matrix(Gstd, tms, xlim_ms, baseWin_this, baseLen_ms, sigWin_N100, ...
            sprintf('%s_%s | Global Standard | ERP (Fz, N100)', subName, expTag), COL);

        subplot(2,2,4);
        plot_snr_curve_fast_safe(Gstd, tms, sigWin_N100, 'neg', nBoot, stepN, minN, ...
            plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, Gstd_stat);
        title(sprintf('%s_%s | Global Standard | SNR (N100)', subName, expTag),'Interpreter','none');

        safe_export(fig4, fullfile(dashOut, sprintf('%s_%s_FIG4_N100.png', subName, expTag)), dpi);
        close(fig4);

        %% ===== write plateau Ns into T_LG (LS/LD/GS/GD for this exp) =====
        row = find(strcmpi(string(T_LG.Subject_ID), string(subName)), 1);
        if ~isempty(row)
            N_LS = plateau_from_curve(NsLS, mLS, plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK);
            N_LD = plateau_from_curve(NsLD, mLD, plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK);
            N_GS = plateau_from_curve(NsGS, mGS, plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK);
            N_GD = plateau_from_curve(NsGD, mGD, plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK);

            T_LG = write_plateau_to_TLG(T_LG, row, expTag, N_LS, N_LD, N_GS, N_GD);
        end

        fprintf('  Saved %s_%s figs -> %s\n', subName, expTag, dashOut);
    end
end

set(groot,'DefaultFigureVisible','on');

%% ===== SAVE updated T_LG (backup + overwrite) =====

Sinfo.T_LG = T_LG;
save(infoFile, '-struct', 'Sinfo');   % overwrite with updated T_LG (keeps other vars)
fprintf('\nUpdated T_LG saved to:\n  %s\n', infoFile);

% Write Excel export 
excelFile = fullfile(outpath,'Info_EEG_LG.xlsx');
write_TLG_excel_exact(T_LG, excelFile);
fprintf('Excel export written: %s\n', excelFile);


%% ===================== LOCAL FUNCTIONS =====================

function T = clear_plateau_this_exp(T, subName, expTag)
row = find(strcmpi(string(T.Subject_ID), string(subName)), 1);
if isempty(row), return; end
if strcmpi(expTag,'E1')
    T.Plateau_E1_LS{row} = [];
    T.Plateau_E1_LD{row} = [];
    T.Plateau_E1_GS{row} = [];
    T.Plateau_E1_GD{row} = [];
else
    T.Plateau_E2_LS{row} = [];
    T.Plateau_E2_LD{row} = [];
    T.Plateau_E2_GS{row} = [];
    T.Plateau_E2_GD{row} = [];
end
end

function X = pick_trials(DATA, ch, idx)
if isempty(idx), X = []; return; end
X = squeeze(DATA(ch,:,idx)); % time x trials
end

function T = ensure_plateau_cols(T)
cols = { ...
    'Plateau_E1_LS','Plateau_E1_LD','Plateau_E1_GS','Plateau_E1_GD', ...
    'Plateau_E2_LS','Plateau_E2_LD','Plateau_E2_GS','Plateau_E2_GD'  ...
    };
for i = 1:numel(cols)
    if ~ismember(cols{i}, T.Properties.VariableNames)
        T.(cols{i}) = repmat({[]}, height(T), 1);
    end
end
end

function T = write_plateau_to_TLG(T, row, expTag, N_LS, N_LD, N_GS, N_GD)
if strcmpi(expTag,'E1')
    T.Plateau_E1_LS{row} = nan_to_empty(N_LS);
    T.Plateau_E1_LD{row} = nan_to_empty(N_LD);
    T.Plateau_E1_GS{row} = nan_to_empty(N_GS);
    T.Plateau_E1_GD{row} = nan_to_empty(N_GD);
else
    T.Plateau_E2_LS{row} = nan_to_empty(N_LS);
    T.Plateau_E2_LD{row} = nan_to_empty(N_LD);
    T.Plateau_E2_GS{row} = nan_to_empty(N_GS);
    T.Plateau_E2_GD{row} = nan_to_empty(N_GD);
end
end

function v = nan_to_empty(x)
if isempty(x) || isnan(x), v = []; else, v = x; end
end

function safe_export(fig, filepath, dpi)
try
    drawnow;
    exportgraphics(fig, filepath, 'Resolution', dpi);
catch
    try, saveas(fig, filepath); end %#ok<TRYNC>
end
end

function EEG = load_or_merge_LG_set(subFold, subName, expTag)
EEG = [];
allSet = dir(fullfile(subFold, sprintf('*auditory_all*LG*%s*%s*.set', subName, expTag)));
if isempty(allSet)
    allSet = dir(fullfile(subFold, sprintf('*all*LG*%s*%s*.set', subName, expTag)));
end
if ~isempty(allSet)
    EEG = pop_loadset('filename', allSet(1).name, 'filepath', allSet(1).folder);
    return;
end

parts = dir(fullfile(subFold, sprintf('*LG*%s*%s*.set', subName, expTag)));
keep = false(size(parts));
for i=1:numel(parts)
    fn = lower(parts(i).name);
    keep(i) = contains(fn,'auditory') && contains(fn,'lg');
end
parts = parts(keep);
if isempty(parts), return; end

EEG = pop_loadset('filename', parts(1).name, 'filepath', parts(1).folder);
for i=2:numel(parts)
    EEG2 = pop_loadset('filename', parts(i).name, 'filepath', parts(i).folder);
    EEG  = pop_mergeset(EEG, EEG2, 1);
end
end

function ch = pick_chan(EEG, want)
ch = find(strcmpi({EEG.chanlocs.labels}, want), 1);
if ~isempty(ch), return; end
if strcmpi(want,'Fz')
    alts = {'FCz','Cz'};
elseif strcmpi(want,'Pz')
    alts = {'CPz','Cz'};
else
    alts = {};
end
ch = [];
for a=1:numel(alts)
    ch = find(strcmpi({EEG.chanlocs.labels}, alts{a}), 1);
    if ~isempty(ch), return; end
end
end

function [locLab, globLab] = get_epoch_local_global_labels(EEG)
nE = EEG.trials;
locLab  = repmat({''}, 1, nE);
globLab = repmat({''}, 1, nE);
for e = 1:nE
    if isfield(EEG.epoch, 'eventlg_sound')
        locLab{e} = pick_first_text(EEG.epoch(e).eventlg_sound);
    end
    if isfield(EEG.epoch, 'eventlg_effect')
        globLab{e} = pick_first_text(EEG.epoch(e).eventlg_effect);
    end
end
end

function s = pick_first_text(v)
s = '';
if isempty(v), return; end
if iscell(v)
    for i=1:numel(v)
        tmp = pick_first_text(v{i});
        if ~isempty(tmp), s = tmp; return; end
    end
    return;
end
if isstring(v), v = char(v); end
if ischar(v), s = strtrim(v); return; end
if isnumeric(v) && isscalar(v), s = num2str(v); end
end

function [bw, x1_ms, baseLen] = get_basewin_from_TLG(T_LG, subName, expTag)
bw = [NaN NaN]; x1_ms = NaN; baseLen = NaN;
if ~istable(T_LG) || ~ismember('Subject_ID', T_LG.Properties.VariableNames), return; end
row = find(strcmpi(string(T_LG.Subject_ID), string(subName)), 1);
if isempty(row), return; end

if strcmpi(expTag,'E1')
    if ~ismember('BaseWin_E1', T_LG.Properties.VariableNames), return; end
    tmp = T_LG.BaseWin_E1{row};
else
    if ~ismember('BaseWin_E2', T_LG.Properties.VariableNames), return; end
    tmp = T_LG.BaseWin_E2{row};
end
if isempty(tmp) || numel(tmp) ~= 2, return; end
tmp = double(tmp(:))';
bw = tmp;
x1_ms = bw(2);
baseLen = bw(2) - bw(1);
end

function [sigWin, peakLat] = pick_peak_window(erp, tms, searchWin, halfWidth, polarity)
idx = tms >= searchWin(1) & tms <= searchWin(2);
if ~any(idx)
    peakLat = mean(searchWin);
    sigWin  = [peakLat-halfWidth, peakLat+halfWidth];
    return;
end
segT = tms(idx);
segY = erp(idx);
if strcmpi(polarity,'neg'), [~,ii] = min(segY); else, [~,ii] = max(segY); end
peakLat = segT(ii);
sigWin  = [peakLat-halfWidth, peakLat+halfWidth];
end

%% ===== BUILDUP plot (same as your group style) =====
function plot_buildup_from_matrix(data, tms, xlim_ms, baseWin, baseLen_ms, sigWin, titleStr, COL)
% plot_buildup_from_matrix (robust + shaded BL/ERP)
% COL.baseRGB / COL.sigRGB / COL.alphaVal control shading colors

if nargin < 8 || isempty(COL)
    COL.baseRGB  = [0.85 1.00 0.85];
    COL.sigRGB   = [0.85 0.92 1.00];
    COL.alphaVal = 0.30;
end

if nargin < 4, baseWin = []; end
if nargin < 5, baseLen_ms = NaN; end
if nargin < 6, sigWin = []; end
if nargin < 7, titleStr = ''; end

tms = tms(:);

baseWin = clip_win(force_row2(unwrap_win(baseWin)), xlim_ms);
sigWin  = clip_win(force_row2(unwrap_win(sigWin)),  xlim_ms);

if isempty(data)
    axis off; text(0.5,0.5,'No trials','HorizontalAlignment','center'); return;
end

nTr = size(data,2);
check = unique(round(linspace(5, nTr, 8)));
check(check < 1) = [];
check(check > nTr) = nTr;
if isempty(check), check = nTr; end

ord = randperm(nTr);

hold on;
ls   = {'-','--',':','-.','-','--',':','-.'};
cmap = lines(max(8, numel(check)));

h = gobjects(numel(check),1);
labels = cell(numel(check),1);
colorIdx = 1;

for i = 1:numel(check)
    N = check(i);
    erp = mean(data(:, ord(1:N)), 2, 'omitnan');

    if N <= 27
        h(i) = plot(tms, erp, 'k', 'LineWidth', 1.2, 'LineStyle', ls{mod(i-1,numel(ls))+1});
    else
        h(i) = plot(tms, erp, 'LineWidth', 2.0, 'Color', cmap(colorIdx,:));
        colorIdx = min(colorIdx+1, size(cmap,1));
    end
    labels{i} = sprintf('n=%d', N);
end

xlim(xlim_ms);
xlabel('Time (ms)'); ylabel('\muV');
grid on; box off;
if ~isempty(titleStr), title(titleStr,'Interpreter','none'); end

% --- Y limits (set to [] to auto) ---
yLimFixed = [-8 8];      % <-- change here (microvolts)

if ~isempty(yLimFixed)
    ylim(yLimFixed);
end

yl = ylim;
yr = yl(2) - yl(1);
if yr == 0, yr = 1; end


% Shaded baseline window
if ~isempty(baseWin)
    xb = [baseWin(1) baseWin(2) baseWin(2) baseWin(1)]';
    yb = [yl(1) yl(1) yl(2) yl(2)]';
    p1 = patch(xb, yb, COL.baseRGB, 'EdgeColor','none', 'FaceAlpha', COL.alphaVal);
    uistack(p1,'bottom');
end

% Shaded ERP window
if ~isempty(sigWin)
    xs = [sigWin(1) sigWin(2) sigWin(2) sigWin(1)]';
    ys = [yl(1) yl(1) yl(2) yl(2)]';
    p2 = patch(xs, ys, COL.sigRGB, 'EdgeColor','none', 'FaceAlpha', COL.alphaVal);
    uistack(p2,'bottom');
end

% BL / ERP text
if ~isempty(baseWin)
    if ~isfinite(baseLen_ms), baseLen_ms = baseWin(2)-baseWin(1); end
    text(mean(baseWin), yl(1) + 0.06*yr, sprintf('BL = %.0f ms', baseLen_ms), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end
if ~isempty(sigWin)
    text(mean(sigWin), yl(1) + 0.06*yr, sprintf('ERP = [%.0f %.0f] ms', sigWin(1), sigWin(2)), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
end

legend(h, labels, 'Location','eastoutside');
end


%% ===== EFFECT BUILDUP (safe) =====
function plot_buildup_effect__SAFE(stdDat, devDat, tms, xlim_ms, baseWin, baseLen_ms, sigWin, titleStr, COL)
% plot_buildup_effect__SAFE
% Difference buildup: mean(dev N) - mean(std N), matched N per side
% Uses COL.baseRGB / COL.sigRGB / COL.alphaVal for shading colors.

% defaults
if nargin < 5, baseWin = []; end
if nargin < 6, baseLen_ms = NaN; end
if nargin < 7, sigWin = []; end
if nargin < 8, titleStr = ''; end

% color defaults if COL not passed
if nargin < 9 || isempty(COL)
    COL.baseRGB  = [0.85 1.00 0.85];
    COL.sigRGB   = [0.85 0.92 1.00];
    COL.alphaVal = 0.30;
end

tms = tms(:);
baseWin = clip_win(force_row2(unwrap_win(baseWin)), xlim_ms);
sigWin  = clip_win(force_row2(unwrap_win(sigWin)),  xlim_ms);

if isempty(stdDat) || isempty(devDat)
    axis off; text(0.5,0.5,'No trials','HorizontalAlignment','center'); return;
end

nS = size(stdDat,2);
nD = size(devDat,2);
nMax = min(nS,nD);
if nMax < 5
    axis off; text(0.5,0.5,'Too few matched trials','HorizontalAlignment','center'); return;
end

check = unique(round(linspace(5, nMax, 8)));
check(check<1)=[]; check(check>nMax)=nMax;
if isempty(check), check = nMax; end

ordS = randperm(nS);
ordD = randperm(nD);

hold on;
ls   = {'-','--',':','-.','-','--',':','-.'};
cmap = lines(max(8,numel(check)));
colorIdx = 1;

h = gobjects(numel(check),1);
labels = cell(numel(check),1);

for i=1:numel(check)
    N = check(i);
    erpS = mean(stdDat(:,ordS(1:N)),2,'omitnan');
    erpD = mean(devDat(:,ordD(1:N)),2,'omitnan');
    diffERP = erpD - erpS;

    if N <= 27
        h(i) = plot(tms, diffERP, 'k', 'LineWidth', 1.2, ...
            'LineStyle', ls{mod(i-1,numel(ls))+1});
    else
        h(i) = plot(tms, diffERP, 'LineWidth', 2.0, 'Color', cmap(colorIdx,:));
        colorIdx = min(colorIdx+1, size(cmap,1));
    end
    labels{i} = sprintf('n=%d/side',N);
end

xlim(xlim_ms);
xlabel('Time (ms)'); ylabel('\muV');
grid on; box off;
if ~isempty(titleStr), title(titleStr,'Interpreter','none'); end


% --- Y limits (set to [] to auto) ---
yLimFixed = [-8 8];   % <-- change here (microvolts)

if ~isempty(yLimFixed)
    ylim(yLimFixed);
end

yl = ylim;
yr = yl(2) - yl(1);
if yr == 0, yr = 1; end


% shaded baseline window
if ~isempty(baseWin)
    patch([baseWin(1) baseWin(2) baseWin(2) baseWin(1)]', ...
          [yl(1) yl(1) yl(2) yl(2)]', ...
          COL.baseRGB, 'EdgeColor','none', 'FaceAlpha', COL.alphaVal);
end

% shaded ERP window
if ~isempty(sigWin)
    patch([sigWin(1) sigWin(2) sigWin(2) sigWin(1)]', ...
          [yl(1) yl(1) yl(2) yl(2)]', ...
          COL.sigRGB, 'EdgeColor','none', 'FaceAlpha', COL.alphaVal);
end

% labels
if ~isempty(baseWin)
    if ~isfinite(baseLen_ms), baseLen_ms = baseWin(2)-baseWin(1); end
    text(mean(baseWin), yl(1)+0.06*yr, sprintf('BL = %.0f ms', baseLen_ms), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end

if ~isempty(sigWin) && all(isfinite(sigWin))
    text(mean(sigWin), yl(1)+0.06*yr, sprintf('ERP = [%.0f %.0f] ms', sigWin(1), sigWin(2)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end

legend(h, labels, 'Location','eastoutside');
end

%% ===== FAST baseline stats (identical noise as concatenation) =====
function S = precompute_base_stats_fixedwin(X, tms, baseWin)
tms = tms(:);
bIdx = tms >= baseWin(1) & tms <= baseWin(2);
if ~any(bIdx)
    error('Baseline window has no samples inside tms. baseWin=[%.1f %.1f].', baseWin(1), baseWin(2));
end
B = X(bIdx, :);
S.nSamp  = sum(bIdx) * ones(1, size(X,2));
S.sum    = sum(B, 1, 'omitnan');
S.sumsq  = sum(B.^2, 1, 'omitnan');
end

function sd = pooled_sd_from_stats(S, idx)
nTot = sum(S.nSamp(idx));
if nTot <= 1, sd = NaN; return; end
sumTot   = sum(S.sum(idx));
sumsqTot = sum(S.sumsq(idx));
var = (sumsqTot - (sumTot.^2)/nTot) / (nTot - 1);
sd  = sqrt(max(var, 0));
end

%% ===== SNR plotters (SAFE) =====
function [Ns, mu] = plot_snr_curve_fast_safe(data, tms, sigWin, pol, nBoot, stepN, minN, ...
    plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, baseStat)

Ns = []; mu = [];

if isempty(data) || size(data,2) < minN || isempty(baseStat)
    axis off;
    nTr = 0; if ~isempty(data), nTr = size(data,2); end
    text(0.5,0.5,sprintf('Not enough trials (n=%d, minN=%d)', nTr, minN), ...
        'HorizontalAlignment','center','Interpreter','none');
    return;
end

if any(isnan(sigWin)) || numel(sigWin)~=2
    axis off;
    text(0.5,0.5,'Signal window undefined','HorizontalAlignment','center','Interpreter','none');
    return;
end

[Ns, mu, lo, hi] = bootstrap_snr_curve_fast(data, tms, sigWin, pol, nBoot, stepN, minN, baseStat);

plot(Ns, mu, 'LineWidth', 2); hold on;
plot(Ns, lo, '--', 'LineWidth', 1);
plot(Ns, hi, '--', 'LineWidth', 1);

yTop = max([hi(:); mu(:)]);
if ~isfinite(yTop), yTop = 1; end

Nplat = compute_plateau(Ns, mu, plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK);
if ~isnan(Nplat)
    xline(Nplat, ':', 'LineWidth', 1.5);
    text(Nplat, yTop, plateau_label(Nplat, plateauMethod, plateauTolBand, slopeTol), ...
        'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9,'Interpreter','none');
else
    text(Ns(1), yTop, sprintf('No plateau (%s)', plateauMethod), ...
        'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9,'Interpreter','none');
end

grid on; box off;
xlabel('# trials'); ylabel('SNR = |signal| / baselineSD');
end

function plot_snr_curve_diff_fast_safe(stdDat, devDat, tms, sigWin, pol, nBoot, stepN, minN, ...
    plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, statStd, statDev)

if isempty(stdDat) || isempty(devDat) || isempty(statStd) || isempty(statDev)
    axis off; text(0.5,0.5,'Not enough trials','HorizontalAlignment','center','Interpreter','none'); return;
end
nMax = min(size(stdDat,2), size(devDat,2));
if nMax < minN
    axis off; text(0.5,0.5,sprintf('Not enough matched trials (n=%d, minN=%d)', nMax, minN), ...
        'HorizontalAlignment','center','Interpreter','none'); return;
end
if any(isnan(sigWin)) || numel(sigWin)~=2
    axis off; text(0.5,0.5,'Signal window undefined','HorizontalAlignment','center','Interpreter','none'); return;
end

[Ns, mu, lo, hi] = bootstrap_snr_curve_diff_fast(stdDat, devDat, tms, sigWin, pol, nBoot, stepN, minN, statStd, statDev);

plot(Ns, mu, 'LineWidth', 2); hold on;
plot(Ns, lo, '--', 'LineWidth', 1);
plot(Ns, hi, '--', 'LineWidth', 1);

yTop = max([hi(:); mu(:)]);
if ~isfinite(yTop), yTop = 1; end

Nplat = compute_plateau(Ns, mu, plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK);
if ~isnan(Nplat)
    xline(Nplat, ':', 'LineWidth', 1.5);
    text(Nplat, yTop, plateau_label(Nplat, plateauMethod, plateauTolBand, slopeTol), ...
        'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9,'Interpreter','none');
else
    text(Ns(1), yTop, sprintf('No plateau (%s)', plateauMethod), ...
        'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9,'Interpreter','none');
end

grid on; box off;
xlabel('# trials per side (matched)'); ylabel('Effect SNR');
end

function [Ns, snrMean, snrLo, snrHi] = bootstrap_snr_curve_fast(data, tms, sigWin, polarity, nBoot, stepN, minN, baseStat)
nTr = size(data,2);
Ns = minN:stepN:nTr;
snrBoot = nan(numel(Ns), nBoot);

for iN = 1:numel(Ns)
    N = Ns(iN);
    for b = 1:nBoot
        idx = randperm(nTr, N);
        erp = mean(data(:,idx), 2, 'omitnan');

        noise  = pooled_sd_from_stats(baseStat, idx);
        signal = signal_from_erp(erp, tms, sigWin, polarity);

        snrBoot(iN,b) = signal / max(noise, eps);
    end
end

snrMean = mean(snrBoot,2,'omitnan');
snrLo   = prctile(snrBoot,  2.5, 2);
snrHi   = prctile(snrBoot, 97.5, 2);
end

function [Ns, snrMean, snrLo, snrHi] = bootstrap_snr_curve_diff_fast(stdDat, devDat, tms, sigWin, polarity, nBoot, stepN, minN, statStd, statDev)
nS = size(stdDat,2);
nD = size(devDat,2);
nMax = min(nS, nD);

Ns = minN:stepN:nMax;
snrBoot = nan(numel(Ns), nBoot);

for iN = 1:numel(Ns)
    N = Ns(iN);
    for b = 1:nBoot
        iS = randperm(nS, N);
        iD = randperm(nD, N);

        erpS = mean(stdDat(:,iS), 2, 'omitnan');
        erpD = mean(devDat(:,iD), 2, 'omitnan');
        diff = erpD - erpS;

        noiseS = pooled_sd_from_stats(statStd, iS);
        noiseD = pooled_sd_from_stats(statDev, iD);
        noise  = sqrt(noiseS.^2 + noiseD.^2);

        signal = signal_from_erp(diff, tms, sigWin, polarity);
        snrBoot(iN,b) = signal / max(noise, eps);
    end
end

snrMean = mean(snrBoot,2,'omitnan');
snrLo   = prctile(snrBoot,  2.5, 2);
snrHi   = prctile(snrBoot, 97.5, 2);
end

function signal = signal_from_erp(erp, tms, sigWin, polarity)
sIdx = tms >= sigWin(1) & tms <= sigWin(2);
seg  = erp(sIdx);
if isempty(seg)
    signal = abs(mean(erp));
    return;
end
if strcmpi(polarity,'neg')
    signal = abs(min(seg));
elseif strcmpi(polarity,'pos')
    signal = abs(max(seg));
else
    signal = abs(mean(seg));
end
end

function N = plateau_from_curve(Ns, mu, method, tolBand, minPtsBand, finalKBand, slopeTol, minPtsSlope, finalKSlope)
if isempty(Ns) || isempty(mu)
    N = NaN;
else
    N = compute_plateau(Ns, mu, method, tolBand, minPtsBand, finalKBand, slopeTol, minPtsSlope, finalKSlope);
end
end

function s = plateau_label(Nplat, method, tolBand, slopeTol)
if strcmpi(method,'band')
    s = sprintf('Plateau = %d (band = %.0f%%)', Nplat, tolBand*100);
else
    s = sprintf('Plateau = %d (slope = %.1f%% / step)', Nplat, slopeTol*100);
end
end

function Nplat = compute_plateau(Ns, mu, method, tolBand, minPtsBand, finalKBand, slopeTol, minPtsSlope, finalKSlope)
Nplat = NaN;
Ns = Ns(:)'; mu = mu(:)';

ok = ~isnan(Ns) & ~isnan(mu);
Ns = Ns(ok); mu = mu(ok);
if numel(mu) < 3, return; end

if strcmpi(method,'band')
    Nplat = plateau_band(Ns, mu, tolBand, minPtsBand, finalKBand);
else
    Nplat = plateau_slope(Ns, mu, slopeTol, minPtsSlope, finalKSlope);
end
end

function Nplat = plateau_band(Ns, mu, tol, minPts, finalK)
Nplat = NaN;
if numel(mu) < minPts, return; end
K = min(finalK, numel(mu));
target = median(mu(end-K+1:end));
absTol = tol * max(abs(target), eps);
bandLo = target - absTol;
bandHi = target + absTol;

for i = 1:(numel(mu)-(minPts-1))
    tail = mu(i:end);
    if all(tail >= bandLo & tail <= bandHi)
        Nplat = Ns(i);
        return;
    end
end
end

function Nplat = plateau_slope(Ns, mu, slopeTol, minPts, finalK)
Nplat = NaN;
if numel(mu) < minPts, return; end
K = min(finalK, numel(mu));
target = median(mu(end-K+1:end));
thr = slopeTol * max(abs(target), eps);

for i = 1:(numel(mu)-(minPts-1))
    tail = mu(i:end);
    d = abs(diff(tail));
    if isempty(d), continue; end
    if all(d <= thr)
        Nplat = Ns(i);
        return;
    end
end
end

function w = unwrap_win(w)
if isempty(w), w = []; return; end
if iscell(w)
    for k = 1:numel(w)
        if ~isempty(w{k}), w = w{k}; break; end
    end
end
if isstring(w), w = char(w); end
end

function w = force_row2(w)
if isempty(w), w = []; return; end
try
    w = double(w(:))';
catch
    w = [];
    return;
end
if numel(w) ~= 2 || any(~isfinite(w))
    w = [];
end
end

function w = clip_win(w, xlim_ms)
if isempty(w) || isempty(xlim_ms), return; end
w(1) = max(w(1), xlim_ms(1));
w(2) = min(w(2), xlim_ms(2));
if w(2) <= w(1), w = []; end
end


