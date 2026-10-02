clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% Calculate_SNR_group_megasubject_TLG.m
% GROUP-as-one (mega-subject) dashboards:
%   - Pools trials across subjects (per experiment) into big time x trials matrices
%   - Same layout as your single-subject dashboards:
%       Fig1: Local Std/Dev (BUILDUP left, SNR right)  [Fz, MMN window]
%       Fig2: Global Std/Dev (BUILDUP left, SNR right) [Pz, P3b window]
%       Fig3: Effects (Dev-Std): BUILDUP(left) + SNR(right) for each effect
%       Fig4: N100 check (X1-aligned, T_LG-based) [Fz]
%
% IMPORTANT (matches your FINAL single-subject logic):
%   - SNR NOISE baseline uses EACH SUBJECT's baseline window from T_LG.BaseWin_E1/E2.
%   - For pooled-condition SNR:
%       noise = SD of ALL baseline samples pooled across selected trials
%       (implemented equivalently via pooled sum/sumsq/nSamples; faster, identical)
%   - For effect (Dev-Std) SNR:
%       noise = sqrt(noiseDev^2 + noiseStd^2)
%
% Plot annotations:
%   - BUILDUP plots show shaded BL (green) + ERP window (blue)
%   - Labels at bottom: "BL = xxx ms" and "ERP = [a b] ms"
%   - Shading baseline window is DISPLAY ONLY (baseWin_this)
%
% Plateau:
%   - Choose 'band' or 'slope'
%
% Saves:
%   - Figures to outFig
%   - Summary table (group plateau N per curve) to MAT + CSV
%
% Miguel Contreras-Altamirano, 2026

%% ===== USER SETTINGS =====
mainpath = 'L:\Downloads\eeglab2023.0\';
outpath  = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';

xlim_ms      = [-2500 700];       % plot x-range (ERP)
baseWin_this = [-2500 -1700];     % baseline shading window (DISPLAY ONLY)
dpi          = 200;

% ===================== NEW: SKIP SUBJECTS LIST =====================
% Put subject folder names you want to exclude (exact folder names, e.g. "S21")
skipSubjects_all = {"S21","S12","S23","S02","S29"};   % example
% ===================================================================

% ERP windows (relative to tone5=0, like your Bekinschtein plots)
MMN.sigWin = [120 220];   MMN.pol = 'neg';   MMN.chan = 'Fz';
P3B.sigWin = [300 600];   P3B.pol = 'pos';   P3B.chan = 'Pz';

% Bootstrap curve settings
nBoot = 100;
minN  = 10;
stepN = 100;     % group is large; use 100 or 200

% Color plots
COL.baseRGB  = [0.5 1 0.5];       % baseline color (green)
COL.sigRGB   = [1.0 0.5 0.0];     % ERP color (blue)
COL.alphaVal = 0.10;              % transparency

% ---------- Plateau criterion ----------
plateauMethod = 'slope';  % 'band' or 'slope'

% Band plateau params:
plateauTolBand = 0.10;    % +/- 10% band
finalK         = 3;       % target = median(last K points)
plateauMinPts  = 3;       % need >=3 points remaining

% Slope plateau params:
slopeTol     = 0.05;      % 5% of target per step
slopeMinPts  = 3;         % need >=3 points remaining
slopeFinalK  = 3;         % target = median(last K points)

% Output
snrOut = fullfile(outpath,'GRAND_AVG');
outFig = fullfile(snrOut,'SNR_Group_level');
if ~exist(outFig,'dir'), mkdir(outFig); end

% Load timing table with BaseWin_E1/E2
infoFile = fullfile(outpath,'Info_EEG_LG.mat');
assert(exist(infoFile,'file')==2, 'Missing: %s', infoFile);
Sinfo = load(infoFile,'T_LG');
T_LG = Sinfo.T_LG;

% N100 relative to X1 (Tone1)
N100.searchRelMS = [60 140];
N100.halfWidthMS = 20;
N100.pol         = 'neg';

%% ===== EEGLAB =====
if exist('eeglab','file') ~= 2
    addpath(mainpath); % ONLY root
end
try, eeglab nogui; catch, eeglab; end
set(groot,'DefaultFigureVisible','off');

%% ===== SUBJECT FOLDERS =====
subDirs = dir(fullfile(outpath,'S*'));
subDirs = subDirs([subDirs.isdir]);
subDirs = subDirs(~ismember({subDirs.name},{'.','..'}));
assert(~isempty(subDirs), 'No S* subject folders found in outpath.');

% ===================== NEW: Apply skipSubjects_all =====================
% Normalize skip list to lower-case strings (robust to "s21" vs "S21")
skipSet = lower(string(skipSubjects_all(:)));
keepMask = true(size(subDirs));
for i = 1:numel(subDirs)
    keepMask(i) = ~ismember(lower(string(subDirs(i).name)), skipSet);
end
subDirs = subDirs(keepMask);

if isempty(subDirs)
    error('After applying skipSubjects_all, no subject folders remain.');
end
% ======================================================================

%% ===== SUMMARY TABLE =====
SUMT = table();

%% ===== RUN BOTH EXPERIMENTS =====
for expTagCell = {'E1','E2'}
    expTag = expTagCell{1};
    fprintf('\n=== GROUP MEGA-SUBJECT | %s ===\n', expTag);

    % -------- collect subject baseline windows from T_LG --------
    baseStart_all = [];
    baseEnd_all   = [];
    baseLen_all   = [];
    x1_all        = []; % baseEnd == X1

    for sd = 1:numel(subDirs)
        subName = subDirs(sd).name;
        [bw, x1_ms, baseLen] = get_basewin_from_TLG(T_LG, subName, expTag);
        if ~any(isnan([bw(:); x1_ms; baseLen]))
            baseStart_all(end+1) = bw(1); %#ok<AGROW>
            baseEnd_all(end+1)   = bw(2); %#ok<AGROW>
            baseLen_all(end+1)   = baseLen; %#ok<AGROW>
            x1_all(end+1)        = x1_ms; %#ok<AGROW>
        end
    end
    assert(~isempty(x1_all), 'No valid BaseWin found in T_LG for %s', expTag);

    % for display label only (median BL length)
    baseLen_ref = median(baseLen_all);

    % for N100 alignment (median X1)
    x1_ref      = median(x1_all);
    baseWin_N100 = [x1_ref - baseLen_ref, x1_ref];

    fprintf('  [%s] BL_len(median)=%.1f ms | X1_ref=%.1f ms | baseWin_N100=[%.1f %.1f]\n', ...
        expTag, baseLen_ref, x1_ref, baseWin_N100(1), baseWin_N100(2));

    % -------- pooled matrices --------
    tms_ref = [];

    Lstd = []; Ldev = [];
    Gstd = []; Gdev = [];

    % baseline stats per trial (pooled) for fast pooled SD
    statLstd = init_stat();  statLdev = init_stat();
    statGstd = init_stat();  statGdev = init_stat();

    % N100 aligned matrices (Fz)
    LstdN100 = [];
    GstdN100 = [];

    % -------- loop subjects --------
    for sd = 1:numel(subDirs)
        subName = subDirs(sd).name;
        subFold = fullfile(subDirs(sd).folder, subDirs(sd).name);

        EEG = load_or_merge_LG_set(subFold, subName, expTag);
        if isempty(EEG) || EEG.trials <= 1
            fprintf('  skip %s %s (missing/empty)\n', subName, expTag);
            continue;
        end

        if isempty(tms_ref)
            tms_ref = EEG.times(:);
        else
            if numel(EEG.times) ~= numel(tms_ref) || any(abs(EEG.times(:)-tms_ref) > 1e-6)
                fprintf('  skip %s %s (times mismatch)\n', subName, expTag);
                continue;
            end
        end

        [bw, x1_ms, ~] = get_basewin_from_TLG(T_LG, subName, expTag);
        if any(isnan([bw(:); x1_ms]))
            fprintf('  skip %s %s (no T_LG BaseWin)\n', subName, expTag);
            continue;
        end

        [locLab, globLab] = get_epoch_local_global_labels(EEG);
        idxLocalStd  = find(strcmpi(locLab,  'Local Standard'));
        idxLocalDev  = find(strcmpi(locLab,  'Local Deviant'));
        idxGlobalStd = find(strcmpi(globLab, 'Global Standard'));
        idxGlobalDev = find(strcmpi(globLab, 'Global Deviant'));

        chF = pick_chan(EEG, 'Fz');
        chP = pick_chan(EEG, 'Pz');
        if isempty(chF) || isempty(chP)
            fprintf('  skip %s %s (missing Fz/Pz fallbacks)\n', subName, expTag);
            continue;
        end

        % baseline mask for this subject (noise baseline)
        bIdx = (tms_ref >= bw(1) & tms_ref <= bw(2));
        if ~any(bIdx)
            fprintf('  skip %s %s (baseline has no samples in epoch)\n', subName, expTag);
            continue;
        end
        nSampPerTrial = sum(bIdx);

        % ---- append Lstd/Ldev (Fz) + baseline stats ----
        if ~isempty(idxLocalStd)
            tmp = squeeze(EEG.data(chF,:,idxLocalStd)); % time x trials
            Lstd = [Lstd, tmp]; %#ok<AGROW>
            statLstd = append_stat(statLstd, tmp(bIdx,:), nSampPerTrial);
        end

        if ~isempty(idxLocalDev)
            tmp = squeeze(EEG.data(chF,:,idxLocalDev));
            Ldev = [Ldev, tmp]; %#ok<AGROW>
            statLdev = append_stat(statLdev, tmp(bIdx,:), nSampPerTrial);
        end

        % ---- append Gstd/Gdev (Pz) + baseline stats ----
        if ~isempty(idxGlobalStd)
            tmp = squeeze(EEG.data(chP,:,idxGlobalStd));
            Gstd = [Gstd, tmp]; %#ok<AGROW>
            statGstd = append_stat(statGstd, tmp(bIdx,:), nSampPerTrial);
        end

        if ~isempty(idxGlobalDev)
            tmp = squeeze(EEG.data(chP,:,idxGlobalDev));
            Gdev = [Gdev, tmp]; %#ok<AGROW>
            statGdev = append_stat(statGdev, tmp(bIdx,:), nSampPerTrial);
        end

        % ---- N100 alignment (Fz) using subject X1 ----
        shift_ms = x1_ref - x1_ms; % shift waveform so X1 lands at x1_ref

        if ~isempty(idxLocalStd)
            tmp  = squeeze(EEG.data(chF,:,idxLocalStd));
            tmpA = shift_trials(tmp, tms_ref, shift_ms);
            LstdN100 = [LstdN100, tmpA]; %#ok<AGROW>
        end
        if ~isempty(idxGlobalStd)
            tmp  = squeeze(EEG.data(chF,:,idxGlobalStd));
            tmpA = shift_trials(tmp, tms_ref, shift_ms);
            GstdN100 = [GstdN100, tmpA]; %#ok<AGROW>
        end
    end

    assert(~isempty(tms_ref), 'No usable subjects found for %s', expTag);

    fprintf('  Trials pooled | LocalStd=%d LocalDev=%d GlobalStd=%d GlobalDev=%d\n', ...
        size(Lstd,2), size(Ldev,2), size(Gstd,2), size(Gdev,2));

    % ===== N100 window from group aligned ERP =====
    sigWin_N100 = [NaN NaN]; peakLat = NaN;
    if ~isempty(LstdN100)
        erpN100 = mean(LstdN100, 2, 'omitnan');
        searchWin = x1_ref + N100.searchRelMS;
        [sigWin_N100, peakLat] = pick_peak_window(erpN100, tms_ref, searchWin, N100.halfWidthMS, 'neg');
    end
    fprintf('  [%s] N100 peak=%.1f ms | ERP=[%.1f %.1f]\n', expTag, peakLat, sigWin_N100(1), sigWin_N100(2));

    % baseline stats for aligned N100 use fixed baseline in aligned coords
    statLstdN100 = init_stat();
    statGstdN100 = init_stat();
    if ~isempty(LstdN100)
        bIdxN = (tms_ref >= baseWin_N100(1) & tms_ref <= baseWin_N100(2));
        if any(bIdxN)
            statLstdN100 = append_stat(statLstdN100, LstdN100(bIdxN,:), sum(bIdxN));
        end
    end
    if ~isempty(GstdN100)
        bIdxN = (tms_ref >= baseWin_N100(1) & tms_ref <= baseWin_N100(2));
        if any(bIdxN)
            statGstdN100 = append_stat(statGstdN100, GstdN100(bIdxN,:), sum(bIdxN));
        end
    end

    %% ================= FIG 1: LOCAL (MMN) =================
    fig1 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

    subplot(2,2,1);
    plot_buildup_from_matrix(Lstd, tms_ref, xlim_ms, baseWin_this, baseLen_ref, MMN.sigWin, ...
        sprintf('GROUP_%s | Local Standard | ERP (Fz)', expTag), COL);

    subplot(2,2,2);
    [NsLS, mLS, NplatLS] = plot_snr_curve_stats_safe(Lstd, tms_ref, MMN.sigWin, MMN.pol, ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statLstd);
    title(sprintf('GROUP_%s | Local Standard | SNR (MMN)', expTag),'Interpreter','none');

    subplot(2,2,3);
    plot_buildup_from_matrix(Ldev, tms_ref, xlim_ms, baseWin_this, baseLen_ref, MMN.sigWin, ...
        sprintf('GROUP_%s | Local Deviant | ERP (Fz)', expTag), COL);

    subplot(2,2,4);
    [NsLD, mLD, NplatLD] = plot_snr_curve_stats_safe(Ldev, tms_ref, MMN.sigWin, MMN.pol, ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statLdev);
    title(sprintf('GROUP_%s | Local Deviant | SNR (MMN)', expTag),'Interpreter','none');

    safe_export(fig1, fullfile(outFig, sprintf('GROUP_%s_FIG1_LOCAL_MMN.png', expTag)), dpi);
    close(fig1);

    %% ================= FIG 2: GLOBAL (P3b) =================
    fig2 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

    subplot(2,2,1);
    plot_buildup_from_matrix(Gstd, tms_ref, xlim_ms, baseWin_this, baseLen_ref, P3B.sigWin, ...
        sprintf('GROUP_%s | Global Standard | ERP (Pz)', expTag), COL);

    subplot(2,2,2);
    [NsGS, mGS, NplatGS] = plot_snr_curve_stats_safe(Gstd, tms_ref, P3B.sigWin, P3B.pol, ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statGstd);
    title(sprintf('GROUP_%s | Global Standard | SNR (P3b)', expTag),'Interpreter','none');

    subplot(2,2,3);
    plot_buildup_from_matrix(Gdev, tms_ref, xlim_ms, baseWin_this, baseLen_ref, P3B.sigWin, ...
        sprintf('GROUP_%s | Global Deviant | ERP (Pz)', expTag), COL);

    subplot(2,2,4);
    [NsGD, mGD, NplatGD] = plot_snr_curve_stats_safe(Gdev, tms_ref, P3B.sigWin, P3B.pol, ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statGdev);
    title(sprintf('GROUP_%s | Global Deviant | SNR (P3b)', expTag),'Interpreter','none');

    safe_export(fig2, fullfile(outFig, sprintf('GROUP_%s_FIG2_GLOBAL_P3B.png', expTag)), dpi);
    close(fig2);

    %% ================= FIG 3: EFFECTS =================
    fig3 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

    % Local effect
    subplot(2,2,1);
    plot_buildup_effect__SAFE(Lstd, Ldev, tms_ref, xlim_ms, baseWin_this, baseLen_ref, MMN.sigWin, ...
        sprintf('GROUP_%s | LOCAL EFFECT (Dev-Std) | ERP (Fz)', expTag), COL);

    subplot(2,2,2);
    [NsLE, mLE, NplatLE] = plot_snr_curve_diff_stats_safe(Lstd, Ldev, tms_ref, MMN.sigWin, MMN.pol, ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statLstd, statLdev);
    title(sprintf('GROUP_%s | LOCAL EFFECT (Dev-Std) | SNR (MMN)', expTag),'Interpreter','none');

    % Global effect
    subplot(2,2,3);
    plot_buildup_effect__SAFE(Gstd, Gdev, tms_ref, xlim_ms, baseWin_this, baseLen_ref, P3B.sigWin, ...
        sprintf('GROUP_%s | GLOBAL EFFECT (Dev-Std) | ERP (Pz)', expTag), COL);

    subplot(2,2,4);
    [NsGE, mGE, NplatGE] = plot_snr_curve_diff_stats_safe(Gstd, Gdev, tms_ref, P3B.sigWin, P3B.pol, ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statGstd, statGdev);
    title(sprintf('GROUP_%s | GLOBAL EFFECT (Dev-Std) | SNR (P3b)', expTag),'Interpreter','none');

    safe_export(fig3, fullfile(outFig, sprintf('GROUP_%s_FIG3_EFFECTS.png', expTag)), dpi);
    close(fig3);

    %% ================= FIG 4: N100 =================
    fig4 = figure('Color','w','Units','pixels','Position',[50 50 1500 900],'Visible','off');

    subplot(2,2,1);
    plot_buildup_from_matrix(LstdN100, tms_ref, xlim_ms, baseWin_this, baseLen_ref, sigWin_N100, ...
        sprintf('GROUP_%s | Local Standard | ERP (Fz, X1-aligned)', expTag), COL);

    subplot(2,2,2);
    [NsN1A, mN1A, NplatN1A] = plot_snr_curve_stats_safe(LstdN100, tms_ref, sigWin_N100, 'neg', ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statLstdN100);
    title(sprintf('GROUP_%s | Local Standard | SNR (N100, aligned)', expTag),'Interpreter','none');

    subplot(2,2,3);
    plot_buildup_from_matrix(GstdN100, tms_ref, xlim_ms, baseWin_this, baseLen_ref, sigWin_N100, ...
        sprintf('GROUP_%s | Global Standard | ERP (Fz, X1-aligned)', expTag), COL);

    subplot(2,2,4);
    [NsN1B, mN1B, NplatN1B] = plot_snr_curve_stats_safe(GstdN100, tms_ref, sigWin_N100, 'neg', ...
        nBoot, stepN, minN, plateauMethod, plateauTolBand, plateauMinPts, finalK, ...
        slopeTol, slopeMinPts, slopeFinalK, statGstdN100);
    title(sprintf('GROUP_%s | Global Standard | SNR (N100, aligned)', expTag),'Interpreter','none');

    safe_export(fig4, fullfile(outFig, sprintf('GROUP_%s_FIG4_N100.png', expTag)), dpi);
    close(fig4);

    fprintf('Saved GROUP figures for %s -> %s\n', expTag, outFig);

    %% ===== store group plateau numbers into summary table =====
    SUMT = [SUMT; make_sum_row(expTag,'Fig1','LocalStd','MMN',  NplatLS, numel_trials(Lstd))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig1','LocalDev','MMN',  NplatLD, numel_trials(Ldev))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig2','GlobalStd','P3b', NplatGS, numel_trials(Gstd))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig2','GlobalDev','P3b', NplatGD, numel_trials(Gdev))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig3','LocalEffect','MMN',  NplatLE, min(numel_trials(Lstd),numel_trials(Ldev)))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig3','GlobalEffect','P3b', NplatGE, min(numel_trials(Gstd),numel_trials(Gdev)))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig4','LocalStd','N100', NplatN1A, numel_trials(LstdN100))]; %#ok<AGROW>
    SUMT = [SUMT; make_sum_row(expTag,'Fig4','GlobalStd','N100', NplatN1B, numel_trials(GstdN100))]; %#ok<AGROW>

end

set(groot,'DefaultFigureVisible','on');

% ---- save summary table ----
sumFileMat = fullfile(outFig, sprintf('GROUP_PLATEAU_SUMMARY_%s.mat', plateauMethod));
save(sumFileMat, 'SUMT');

try
    sumFileCsv = fullfile(outFig, sprintf('GROUP_PLATEAU_SUMMARY_%s.csv', plateauMethod));
    writetable(SUMT, sumFileCsv);
catch
end

fprintf('\nSaved group summary table:\n  %s\n', sumFileMat);

%% ========================= LOCAL HELPERS =========================

function n = numel_trials(X)
if isempty(X), n = 0; else, n = size(X,2); end
end

function rowT = make_sum_row(expTag, figName, condition, comp, Nplat, Ntrials)
if isempty(Nplat) || isnan(Nplat), NplatStr = "[]"; else, NplatStr = string(Nplat); end
rowT = table( ...
    string(expTag), string(figName), string(condition), string(comp), ...
    Ntrials, NplatStr, ...
    'VariableNames', {'Experiment','Figure','Condition','Component','Ntrials','PlateauN'} );
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

function Xs = shift_trials(X, tms, shift_ms)
if isempty(X), Xs = X; return; end
tms = tms(:);
Xs = nan(size(X));
tQuery = tms - shift_ms;
for tr = 1:size(X,2)
    y = X(:,tr);
    Xs(:,tr) = interp1(tms, y, tQuery, 'linear', NaN);
end
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

%% ===== baseline stats (fast pooled SD) =====
function S = init_stat()
S.sum = []; S.sumsq = []; S.nSamp = [];
end

function S = append_stat(S, B, nSampPerTrial)
% B: baseline samples (nBaseSamp x nTrialsBlock)
if isempty(B), return; end
S.sum   = [S.sum,   sum(B,1,'omitnan')];
S.sumsq = [S.sumsq, sum(B.^2,1,'omitnan')];
S.nSamp = [S.nSamp, nSampPerTrial * ones(1,size(B,2))];
end

function sd = pooled_sd_from_stats(S, idx)
% idx indexes trials in the pooled matrix
if isempty(idx) || isempty(S.sum), sd = NaN; return; end
nTot = sum(S.nSamp(idx));
if nTot <= 1, sd = NaN; return; end
sumTot   = sum(S.sum(idx));
sumsqTot = sum(S.sumsq(idx));
var = (sumsqTot - (sumTot.^2)/nTot) / (nTot - 1);
sd  = sqrt(max(var, 0));
end

%% ===== SNR curves (stats-safe) =====
function [Ns, mu, Nplat] = plot_snr_curve_stats_safe(data, tms, sigWin, pol, nBoot, stepN, minN, ...
    plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, stat)

Ns=[]; mu=[]; Nplat=NaN;

if isempty(data) || size(data,2) < minN || isempty(stat.sum)
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

[Ns, mu, lo, hi] = bootstrap_snr_curve_stats(data, tms, sigWin, pol, nBoot, stepN, minN, stat);

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

function [Ns, mu, Nplat] = plot_snr_curve_diff_stats_safe(stdDat, devDat, tms, sigWin, pol, nBoot, stepN, minN, ...
    plateauMethod, plateauTolBand, plateauMinPts, finalK, slopeTol, slopeMinPts, slopeFinalK, statStd, statDev)

Ns=[]; mu=[]; Nplat=NaN;

if isempty(stdDat) || isempty(devDat)
    axis off; text(0.5,0.5,'No trials','HorizontalAlignment','center','Interpreter','none'); return;
end
nMax = min(size(stdDat,2), size(devDat,2));
if nMax < minN || isempty(statStd.sum) || isempty(statDev.sum)
    axis off; text(0.5,0.5,sprintf('Not enough matched trials (n=%d, minN=%d)', nMax, minN), ...
        'HorizontalAlignment','center','Interpreter','none'); return;
end
if any(isnan(sigWin)) || numel(sigWin)~=2
    axis off; text(0.5,0.5,'Signal window undefined','HorizontalAlignment','center','Interpreter','none'); return;
end

[Ns, mu, lo, hi] = bootstrap_snr_curve_diff_stats(stdDat, devDat, tms, sigWin, pol, nBoot, stepN, minN, statStd, statDev);

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

function [Ns, mu, lo, hi] = bootstrap_snr_curve_stats(data, tms, sigWin, polarity, nBoot, stepN, minN, stat)
nTr = size(data,2);
Ns = minN:stepN:nTr;
snrBoot = nan(numel(Ns), nBoot);

for iN = 1:numel(Ns)
    N = Ns(iN);
    for b = 1:nBoot
        idx = randperm(nTr, N);
        erp = mean(data(:,idx), 2, 'omitnan');

        noise  = pooled_sd_from_stats(stat, idx);
        signal = signal_from_erp(erp, tms, sigWin, polarity);

        snrBoot(iN,b) = signal / max(noise, eps);
    end
end

mu = mean(snrBoot,2,'omitnan');
lo = prctile(snrBoot,  2.5, 2);
hi = prctile(snrBoot, 97.5, 2);
end

function [Ns, mu, lo, hi] = bootstrap_snr_curve_diff_stats(stdDat, devDat, tms, sigWin, polarity, nBoot, stepN, minN, statStd, statDev)
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

        erpS = mean(stdDat(:,iS),2,'omitnan');
        erpD = mean(devDat(:,iD),2,'omitnan');
        diff = erpD - erpS;

        noiseS = pooled_sd_from_stats(statStd, iS);
        noiseD = pooled_sd_from_stats(statDev, iD);
        noise  = sqrt(noiseS.^2 + noiseD.^2);

        signal = signal_from_erp(diff, tms, sigWin, polarity);
        snrBoot(iN,b) = signal / max(noise, eps);
    end
end

mu = mean(snrBoot,2,'omitnan');
lo = prctile(snrBoot,  2.5, 2);
hi = prctile(snrBoot, 97.5, 2);
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

%% ===== plateau =====
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

%% ===== BUILDUP plots (same style as your stable versions) =====
function plot_buildup_from_matrix(data, tms, xlim_ms, baseWin, baseLen_ms, sigWin, titleStr, COL)
if nargin < 4, baseWin = []; end
if nargin < 5, baseLen_ms = NaN; end
if nargin < 6, sigWin = []; end
if nargin < 7, titleStr = ''; end

if nargin < 8 || isempty(COL)
    COL.baseRGB  = [0.5 1 0.5];       % baseline color (green)
    COL.sigRGB   = [1.0 0.5 0.0];     % ERP color (blue)
    COL.alphaVal = 0.10;              % transparency
end

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
yLimFixed = [-8 8];   % <-- change here (microvolts)

if ~isempty(yLimFixed)
    ylim(yLimFixed);
end

yl = ylim;
yr = yl(2) - yl(1);
if yr == 0, yr = 1; end

if ~isempty(baseWin)
    patch([baseWin(1) baseWin(2) baseWin(2) baseWin(1)]', [yl(1) yl(1) yl(2) yl(2)]', ...
        COL.baseRGB,'EdgeColor','none','FaceAlpha',COL.alphaVal);
end
if ~isempty(sigWin)
    patch([sigWin(1) sigWin(2) sigWin(2) sigWin(1)]', [yl(1) yl(1) yl(2) yl(2)]', ...
        COL.sigRGB,'EdgeColor','none','FaceAlpha',COL.alphaVal);
end

if ~isempty(baseWin)
    if ~isfinite(baseLen_ms), baseLen_ms = baseWin(2)-baseWin(1); end
    text(mean(baseWin), yl(1)+0.06*yr, sprintf('BL = %.0f ms', baseLen_ms), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end
if ~isempty(sigWin) && all(isfinite(sigWin))
    text(mean(sigWin), yl(1)+0.06*yr, sprintf('ERP = [%.0f %.0f] ms', sigWin(1), sigWin(2)), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
end

legend(h, labels, 'Location','eastoutside');
end

function plot_buildup_effect__SAFE(stdDat, devDat, tms, xlim_ms, baseWin, baseLen_ms, sigWin, titleStr, COL)
if nargin < 5, baseWin = []; end
if nargin < 6, baseLen_ms = NaN; end
if nargin < 7, sigWin = []; end
if nargin < 8, titleStr = ''; end

tms = tms(:);
baseWin = clip_win(force_row2(unwrap_win(baseWin)), xlim_ms);
sigWin  = clip_win(force_row2(unwrap_win(sigWin)),  xlim_ms);

if isempty(stdDat) || isempty(devDat)
    axis off; text(0.5,0.5,'No trials','HorizontalAlignment','center'); return;
end

nMax = min(size(stdDat,2), size(devDat,2));
if nMax < 5
    axis off; text(0.5,0.5,'Too few matched trials','HorizontalAlignment','center'); return;
end

check = unique(round(linspace(5, nMax, 8)));
check(check<1)=[]; check(check>nMax)=nMax;
if isempty(check), check = nMax; end

ordS = randperm(size(stdDat,2));
ordD = randperm(size(devDat,2));

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
        h(i) = plot(tms, diffERP, 'k', 'LineWidth', 1.2, 'LineStyle', ls{mod(i-1,numel(ls))+1});
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

if ~isempty(baseWin)
    patch([baseWin(1) baseWin(2) baseWin(2) baseWin(1)]', [yl(1) yl(1) yl(2) yl(2)]', ...
        COL.baseRGB,'EdgeColor','none','FaceAlpha',COL.alphaVal);
end
if ~isempty(sigWin) && all(isfinite(sigWin))
    patch([sigWin(1) sigWin(2) sigWin(2) sigWin(1)]', [yl(1) yl(1) yl(2) yl(2)]', ...
        COL.sigRGB,'EdgeColor','none','FaceAlpha',COL.alphaVal);
end

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
try, w = double(w(:))'; catch, w = []; return; end
if numel(w) ~= 2 || any(~isfinite(w)), w = []; end
end

function w = clip_win(w, xlim_ms)
if isempty(w) || isempty(xlim_ms), return; end
w(1) = max(w(1), xlim_ms(1));
w(2) = min(w(2), xlim_ms(2));
if w(2) <= w(1), w = []; end
end