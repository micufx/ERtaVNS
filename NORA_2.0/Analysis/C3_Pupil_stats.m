clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% PUPIL STATS (Mazancieux et al., 2023-style) + GRAND-AVG FIGURES WITH TOP SIGNIFICANCE BARS
%
% PURPOSE
%   - Load preprocessed pupil MAT files for each subject (E1 and E2).
%   - Compute subject-level mean pupil time-courses for:
%       (1) Local effect:  Local Standard (LS) vs Local Deviant (LD)   [from dataTableClean.sound]
%       (2) Global effect: Global Standard (GS) vs Global Deviant (GD) [from Frequent vs Rare mapping of dataTableClean.type]
%       (3) Interaction:   (FLD - FLS) - (RLD - RLS)
%   - Run 1D cluster-based permutation tests (FWE-corrected) across time.
%   - Plot grand-averages (mean ± SEM across subjects) with a TOP significance bar.
%
% CRITICAL NOTE (FIX)
%   Your original grand-average LS/LD curves are computed DIRECTLY from:
%       dataTableClean.sound  ("Local Standard" vs "Local Deviant")
%   A previous approach reconstructing LS/LD as mean(FLS,RLS) / mean(FLD,RLD)
%   would weight subconditions equally and can differ from trial-wise averaging.
%
%   This script computes:
%     - LS and LD subject means DIRECTLY from "sound" (MATCHES YOUR GRAND-AVG).
%     - GS and GD subject means from trial-wise mapping of TYPE (Frequent/Rare).
%     - Interaction from TYPE means: (FLD-FLS) - (RLD-RLS).
%
% OUTPUT FIGURES (per experiment)
%   1) Local (LS vs LD)  + TOP sig bar + "pFWE < 0.05"
%   2) Type (FLS,RLS,FLD,RLD) + TOP sig bar (THIS BAR = Global effect GS vs GD)
%   3) 2-panel (Local + Global) + TOP sig bars
%
% INTERPRETATION OF THE BLUE BAR
%   - Local figure:  bar corresponds ONLY to LS vs LD cluster-permutation test.
%   - Type figure:   bar corresponds ONLY to GS vs GD (Frequent vs Rare) test,
%                    NOT "RLD vs the other three". It is a collapsed global contrast.
%   - Local+Global:  top bar = LS vs LD ; bottom bar = GS vs GD.
%
% INTERACTION
% (FLD - FLS) - (RLD - RLS)
%
% Miguel Contreras Altamirano, 2026

%% ===================== USER SETTINGS =====================
outpath  = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';
useField = 'PUPIL_Z';

saveFigures   = true;
saveStatsMat  = true;

grandOutFolder = fullfile(outpath, 'GRAND_AVG');
if ~exist(grandOutFolder, 'dir'); mkdir(grandOutFolder); end

excludeSubs = ["S21","S12","S23","S02","S29"];   % keep [] for none

% ===================== CLUSTER TEST SETTINGS =====================
nPerm          = 10000;
alphaCluster   = 0.05;    % cluster-forming alpha (two-tailed)
alphaFWE       = 0.05;    % cluster-level alpha (FWE)
rngSeed        = 1;
rng(rngSeed, 'twister');  % reproducible random trial balancing

%% ===================== TESTING WINDOW (NEW) =====================
% Limits (ms) for the CLUSTER PERMUTATION TEST ONLY.
% The plots can still show the full epoch; the stats only use this window.
statsWinMS = [0 2000];

% If true, x-axis will be zoomed to statsWinMS for plots as well.
% If false, plots show the full epoch (recommended).
zoomPlotsToStatsWin = false;

%% ===================== SIGNIFICANCE BAR STYLE (NEW) =====================
% Add some space above the data so the bar doesn't cover the curves.
headroomFrac     = 0.08;   % 10% of y-range added to the TOP
barFracFromTop   = 0.06;   % bar placed 2% below top (after headroom)
pTxtYOffsetFrac  = 0.015;  % p-text slightly above the bar

%% ===================== SSI BIN COLORS (EDIT HERE) =====================
% These control the Bekinschtein-like SSI MATRIX colors.
% Rule:
%   - You have nBins = numel(pBins_matrix)
%   - You MUST provide (nBins+1) colors:
%       row 1 = level 0  ("none")
%       row 2 = level 1  (<= pBins_matrix(1))
%       ...
%       row nBins+1 = level nBins (<= pBins_matrix(end))

SSI_BINS.pBins_matrix = [0.05 0.01 0.005 0.001 0.0005 0.0001];  % EDIT
SSI_BINS.noneLabel    = ">0.05";
SSI_BINS.binLabel     = "p-values";

% ---- EDIT COLORS HERE (RGB rows) ----
% Smooth palette: light "none" + (dark blue -> red) for bins 1..nBins
a=[0 0 0];  % very light (none)
blueDark = [33,113,181]/255;
redDark  = [178,24,43]/255;

nBins = numel(SSI_BINS.pBins_matrix);      % here = 6
gradBR = interp1([0 1], [blueDark; redDark], linspace(0,1,nBins), 'pchip'); % 6 colors

SSI_BINS.binColors = [ ...
    a;        % level 0 = none
    gradBR;   % levels 1..6 = dark blue -> red (smooth)
];

%% ===== BINNED MATRIX EXPORT (optional) =====
doBinnedMatrix = true;

cfgBIN = struct();
cfgBIN.binSizeMS = 100;
cfgBIN.binStepMS = 50;
cfgBIN.winMS     = [0 2000];
cfgBIN.doFWE     = true;   % TRUE or FALSE
cfgBIN.fweMethod = "holm";  % or "bonferroni"
cfgBIN.alphaFWE  = 0.05;

pupil_export_binned_subject_matrix(string(outpath), string(grandOutFolder), cfgBIN, excludeSubs, string(useField));


%% ===================== SIGNIFICANCE BAR TEXT SWITCH (NEW) =====================
% If false: bars are still shown, but NO "pFWE < ..." text is added anywhere.
% In the TYPE compare figure (E1 vs E2), we will still show "Local"/"Global" labels.
pText = false;   % true / false


%% ===================== FIND PARTICIPANTS =====================
if ~exist(outpath,'dir')
    error('outpath does not exist or is not mounted: %s', outpath);
end

dSubs = dir(fullfile(outpath, 'S*'));
dSubs = dSubs([dSubs.isdir]);

names = string({dSubs.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));

keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+', 'once'));
subNames = string(names(keep));

if isempty(subNames)
    error('No participant folders found under: %s  (expected folders like S01, S02, ...)', outpath);
end

if ~isempty(excludeSubs)
    beforeN = numel(subNames);
    subNames = subNames(~ismember(subNames, excludeSubs));
    afterN  = numel(subNames);
    fprintf('Excluded %d participant(s): %s\n', beforeN-afterN, strjoin(excludeSubs, ", "));
end

fprintf('Participants included for analysis: %d\n', numel(subNames));
if isempty(subNames)
    error('All participants were excluded. Nothing to analyze.');
end

%% ===================== CONDITIONS (MUST MATCH YOUR LABELS) =====================
soundOrder  = ["Local Standard", "Local Deviant"];

typeOrder  = [ ...
    "Frequent Local Standard", ...
    "Rare Local Standard", ...
    "Frequent Local Deviant", ...
    "Rare Local Deviant"];

effectOrder = ["Global Standard","Global Deviant"]; % computed from TYPE (Frequent/Rare)

%% ===================== MAIN LOOP: E1 & E2 =====================
STATS = struct();
GRAND = struct();
SSI = struct();

for eidx = 1:2
    thisExp = sprintf('E%d', eidx);

    fprintf('\n====================================================\n');
    fprintf('PUPIL STATS + GRAND FIGS (FIXED LS/LD): %s\n', thisExp);
    fprintf('====================================================\n');

    relTimes_ref = [];
    subjIDs      = strings(0,1);

    % Subject mean matrices (nSub x nTime)
    M_LS  = [];
    M_LD  = [];
    M_GS  = [];
    M_GD  = [];

    M_FLS = [];
    M_RLS = [];
    M_FLD = [];
    M_RLD = [];

    % trial totals across included subjects
    tot_LS = 0; tot_LD = 0;
    tot_GS = 0; tot_GD = 0;

    tot_FLS = 0; tot_RLS = 0; tot_FLD = 0; tot_RLD = 0;


    % Mean index stats
    SSI.(thisExp) = struct();
    SSI.(thisExp).times    = relTimes_ref; % will set after first subject
    SSI.(thisExp).local    = [];
    SSI.(thisExp).global   = [];
    SSI.(thisExp).pLoc     = [];
    SSI.(thisExp).pGlo     = [];
    SSI.(thisExp).subjects = strings(0,1);


    % ===== Bekinschtein-like individual stats (SSI) from TRIALS =====
cfgSSI.minConsec_ind = 20;  
cfgSSI.alpha_ind     = 0.05;  % first criterion p<=0.05
cfgSSI.post0_only    = true;  % match paper figures
cfgSSI.exportDPI     = 300;

cfgSSI.pBins_matrix = SSI_BINS.pBins_matrix;
cfgSSI.binColors    = SSI_BINS.binColors;
cfgSSI.noneLabel    = SSI_BINS.noneLabel;
cfgSSI.binLabel     = SSI_BINS.binLabel;

% ===== INDIVIDUAL TRIAL-LEVEL TEST OPTIONS =====
% Applies ONLY to the dark/red Bekinschtein-like individual SSI matrices.
% It does NOT affect the group-level cluster permutation tests.
cfgSSI.testMethod    = "student";      % "welch" or "student"
cfgSSI.balanceTrials = true;         % true = randomly equalize trial numbers
cfgSSI.balanceMode   = "elementary"; % "elementary" or "contrast"


    %% ---------- SUBJECT LOOP ----------
    for s = 1:numel(subNames)
        participant = subNames(s);
        subFolder   = fullfile(outpath, char(participant));

        matFile = fullfile(subFolder, sprintf('%s_%s_pupil_preproc.mat', participant, thisExp));
        if ~exist(matFile, 'file')
            fprintf('  [%s | %s] MAT not found -> skipping\n', participant, thisExp);
            continue;
        end

        fprintf('\n  [%s | %s] Loading: %s\n', participant, thisExp, matFile);
        load(matFile, 'dataTableClean', 'relTimes');

        if ~exist('dataTableClean','var') || isempty(dataTableClean) || height(dataTableClean)==0
            fprintf('  [%s | %s] dataTableClean empty -> skipping\n', participant, thisExp);
            clear dataTableClean relTimes
            continue;
        end

% Check time vector consistency
if isempty(relTimes_ref)
    relTimes_ref = relTimes(:)';
    fprintf('  [%s | %s] relTimes set as reference (%d samples)\n', participant, thisExp, numel(relTimes_ref));
else
    if numel(relTimes) ~= numel(relTimes_ref) || any(relTimes(:)' ~= relTimes_ref)
        error('Time vector mismatch in %s (%s). Ensure all files use same relTimes.', participant, thisExp);
    end
end

% (Optional) now that relTimes_ref is known, keep SSI times synced
SSI.(thisExp).times = relTimes_ref;

        % ---------- Trial counts ----------
        soundVec = string(dataTableClean.sound);
        nLS = sum(soundVec == soundOrder(1));
        nLD = sum(soundVec == soundOrder(2));

        typeVec = string(dataTableClean.type);
        nFLS = sum(typeVec == typeOrder(1));
        nRLS = sum(typeVec == typeOrder(2));
        nFLD = sum(typeVec == typeOrder(3));
        nRLD = sum(typeVec == typeOrder(4));

        % Trial-wise global effect mapping (Frequent = GS, Rare = GD)
        effectVec = strings(height(dataTableClean),1);
        effectVec(startsWith(typeVec, "Frequent")) = "Global Standard";
        effectVec(startsWith(typeVec, "Rare"))     = "Global Deviant";

        nGS = sum(effectVec == effectOrder(1));
        nGD = sum(effectVec == effectOrder(2));

        fprintf('  [%s | %s] SOUND counts -> LS=%d | LD=%d\n', participant, thisExp, nLS, nLD);
        fprintf('  [%s | %s] TYPE  counts -> FLS=%d | RLS=%d | FLD=%d | RLD=%d\n', participant, thisExp, nFLS, nRLS, nFLD, nRLD);
        fprintf('  [%s | %s] GLOBAL counts -> GS=%d | GD=%d\n', participant, thisExp, nGS, nGD);


        %% Matrix individual stats

        % Local: LS vs LD (from SOUND)
        [ssiLoc, pLoc] = pupil_individual_pmap_bekinschtein( ...
            dataTableClean, useField, relTimes_ref, ...
            string(dataTableClean.sound), "Local Standard", "Local Deviant", cfgSSI);

        % Global: GS vs GD (from effectVec you already created)
        [ssiGlo, pGlo] = pupil_individual_pmap_bekinschtein( ...
            dataTableClean, useField, relTimes_ref, ...
            effectVec, "Global Standard", "Global Deviant", cfgSSI);

        % Store (you need to initialize SSI struct before the subject loop)
        SSI.(thisExp).times = relTimes_ref;
        SSI.(thisExp).local(end+1,:)  = ssiLoc; %#ok<SAGROW>
        SSI.(thisExp).global(end+1,:) = ssiGlo; %#ok<SAGROW>
        SSI.(thisExp).pLoc(end+1,:)   = pLoc;   %#ok<SAGROW>
        SSI.(thisExp).pGlo(end+1,:)   = pGlo;   %#ok<SAGROW>
        SSI.(thisExp).subjects(end+1,1) = participant; %#ok<SAGROW>


        %%

        % ---------- SUBJECT MEANS ----------
        % IMPORTANT FIX: LS/LD come DIRECTLY from "sound" (matches your grand-average)
        mLS = subject_mean_curve(dataTableClean, useField, 'sound', soundOrder(1));
        mLD = subject_mean_curve(dataTableClean, useField, 'sound', soundOrder(2));

        % Global effect means from trial-wise effectVec (NOT mean(FLS,FLD), etc.)
        Ttmp = dataTableClean;
        Ttmp.effect = effectVec;
        mGS = subject_mean_curve(Ttmp, useField, 'effect', effectOrder(1));
        mGD = subject_mean_curve(Ttmp, useField, 'effect', effectOrder(2));

        % Type means (for TYPE plot + interaction)
        mFLS = subject_mean_curve(dataTableClean, useField, 'type', typeOrder(1));
        mRLS = subject_mean_curve(dataTableClean, useField, 'type', typeOrder(2));
        mFLD = subject_mean_curve(dataTableClean, useField, 'type', typeOrder(3));
        mRLD = subject_mean_curve(dataTableClean, useField, 'type', typeOrder(4));

        % If subject contributes nothing at all -> skip
        if all(isnan(mLS)) && all(isnan(mLD)) && all(isnan(mGS)) && all(isnan(mGD)) ...
                && all(isnan(mFLS)) && all(isnan(mRLS)) && all(isnan(mFLD)) && all(isnan(mRLD))
            fprintf('  [%s | %s] No usable trials -> skipping subject.\n', participant, thisExp);
            clear dataTableClean relTimes Ttmp effectVec
            continue;
        end

        % Append matrices
        M_LS(end+1,:) = mLS; %#ok<SAGROW>
        M_LD(end+1,:) = mLD; %#ok<SAGROW>
        M_GS(end+1,:) = mGS; %#ok<SAGROW>
        M_GD(end+1,:) = mGD; %#ok<SAGROW>

        M_FLS(end+1,:) = mFLS; %#ok<SAGROW>
        M_RLS(end+1,:) = mRLS; %#ok<SAGROW>
        M_FLD(end+1,:) = mFLD; %#ok<SAGROW>
        M_RLD(end+1,:) = mRLD; %#ok<SAGROW>

        subjIDs(end+1,1) = participant; %#ok<SAGROW>

        % Update totals (titles)
        tot_LS = tot_LS + nLS; tot_LD = tot_LD + nLD;
        tot_GS = tot_GS + nGS; tot_GD = tot_GD + nGD;

        tot_FLS = tot_FLS + nFLS; tot_RLS = tot_RLS + nRLS;
        tot_FLD = tot_FLD + nFLD; tot_RLD = tot_RLD + nRLD;

        clear dataTableClean relTimes Ttmp effectVec
    end

    if isempty(subjIDs)
        warning('No valid subjects for %s. Skipping.', thisExp);
        continue;
    end

    Nsub = numel(subjIDs);
    fprintf('\n[%s] Subjects included in stats: %d\n', thisExp, Nsub);

    relTimes_ref = relTimes_ref(:)';

    %% ---------- GRAND MEAN + SEM ----------
    % SOUND grand (LS vs LD)
    [mu_sound, sem_sound] = mean_sem_over_subjects({M_LS, M_LD});

    % TYPE grand (FLS,RLS,FLD,RLD)
    [mu_type, sem_type]   = mean_sem_over_subjects({M_FLS, M_RLS, M_FLD, M_RLD});

    % EFFECT grand (GS vs GD)
    [mu_eff, sem_eff]     = mean_sem_over_subjects({M_GS, M_GD});

    Ntr_local  = [tot_LS; tot_LD];
    Ntr_global = [tot_GS; tot_GD];

    %% ---------- CLUSTER PERMUTATION TESTS (WINDOWED) ----------
    rng(rngSeed);

    fprintf('[%s] Running cluster permutation tests (nPerm=%d) in window [%d..%d] ms...\n', ...
        thisExp, nPerm, statsWinMS(1), statsWinMS(2));

    % Local effect (paired): LS vs LD (WINDOWED)
    [maskLocal, clustersLocal] = clusterperm_paired_1d_windowed( ...
        M_LS, M_LD, relTimes_ref, statsWinMS, alphaCluster, alphaFWE, nPerm);

    % Global effect (paired): GS vs GD (WINDOWED)
    [maskGlobal, clustersGlobal] = clusterperm_paired_1d_windowed( ...
        M_GS, M_GD, relTimes_ref, statsWinMS, alphaCluster, alphaFWE, nPerm);

% Interaction (local x global) (WINDOWED)
% Same statistical family as the blue bars:
% paired cluster-permutation across subjects.
%
% We compare the LOCAL effect under the two GLOBAL contexts:
%   A_int = FLD - FLS   (local effect when pattern is Frequent / Global Standard)
%   B_int = RLD - RLS   (local effect when pattern is Rare     / Global Deviant)
%
% This is algebraically equivalent to testing:
%   (FLD - FLS) - (RLD - RLS)  against 0
% but expressed as a paired test, matching the blue-bar logic.

A_int = M_FLD - M_FLS;
B_int = M_RLD - M_RLS;

[maskInt, clustersInt] = clusterperm_paired_1d_windowed( ...
    A_int, B_int, relTimes_ref, statsWinMS, alphaCluster, alphaFWE, nPerm);

    fprintf('[%s] Done. Sig timepoints: Local=%d | Global=%d | Interaction=%d\n', ...
        thisExp, nnz(maskLocal), nnz(maskGlobal), nnz(maskInt));

    %% ---------- STORE OUTPUTS ----------
    STATS.(thisExp).subjects      = subjIDs;
    STATS.(thisExp).relTimes      = relTimes_ref;
    STATS.(thisExp).nPerm         = nPerm;
    STATS.(thisExp).alphaCluster  = alphaCluster;
    STATS.(thisExp).alphaFWE      = alphaFWE;
    STATS.(thisExp).statsWinMS    = statsWinMS;

    STATS.(thisExp).maskLocal     = maskLocal;
    STATS.(thisExp).maskGlobal    = maskGlobal;
    STATS.(thisExp).maskInteract  = maskInt;

    STATS.(thisExp).clustersLocal    = clustersLocal;
    STATS.(thisExp).clustersGlobal   = clustersGlobal;
    STATS.(thisExp).clustersInteract = clustersInt;

    GRAND.(thisExp).relTimes      = relTimes_ref;
    GRAND.(thisExp).mu_sound      = mu_sound;
    GRAND.(thisExp).sem_sound     = sem_sound;
    GRAND.(thisExp).mu_type       = mu_type;
    GRAND.(thisExp).sem_type      = sem_type;
    GRAND.(thisExp).mu_effect     = mu_eff;
    GRAND.(thisExp).sem_effect    = sem_eff;

    GRAND.(thisExp).Nsub          = Nsub;
    GRAND.(thisExp).Ntr_local     = Ntr_local;
    GRAND.(thisExp).Ntr_global    = Ntr_global;


    if saveStatsMat
        outMat = fullfile(grandOutFolder, sprintf('STATS_%s_pupil_clusterperm.mat', thisExp));
        STATS_this = STATS.(thisExp); %#ok<NASGU>
        GRAND_this = GRAND.(thisExp); %#ok<NASGU>
        save(outMat,'STATS_this','GRAND_this','-v7.3');
        fprintf('[%s] Saved stats MAT: %s\n', thisExp, outMat);
    end


    % ---------- INTERACTION SUMMARY FOR PAPER ----------
STATS.(thisExp).interactionSummary = summarize_interaction_for_paper( ...
    A_int, B_int, relTimes_ref, statsWinMS, clustersInt, alphaFWE, thisExp);

fprintf('[%s] %s\n', thisExp, STATS.(thisExp).interactionSummary.sentence);



    %% ===================== PLOTS (WITH TOP SIG BARS + TEXT) =====================

pTxtLocal  = ptext_simple(clustersLocal, alphaFWE);
pTxtGlobal = ptext_simple(clustersGlobal, alphaFWE);

% ---- NEW: global switch to show/hide p-text everywhere ----
if ~pText
    pTxtLocal  = "";
    pTxtGlobal = "";
end

    % --- Figure 1: LOCAL (LS vs LD) + local sig bar ---
fig1 = plot_grand_mean_sem_local_with_sig( ...
    relTimes_ref, soundOrder, mu_sound, sem_sound, ...
    thisExp, Nsub, Ntr_local, maskLocal, maskInt, pTxtLocal, ...
    statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac);

    if saveFigures
        outPng = fullfile(grandOutFolder, sprintf('STATS_%s_meanSEM_bySound_sig.png', thisExp));
        exportgraphics(fig1, outPng, 'BackgroundColor','white', 'Resolution',300);
        close(fig1);
        fprintf('[%s] Saved: %s\n', thisExp, outPng);
    end

    % --- Figure 2: TYPE (4 lines) + global sig bar (GS vs GD) ---
    fig2 = plot_grand_mean_sem_by_type_with_sig( ...
        relTimes_ref, typeOrder, mu_type, sem_type, ...
        thisExp, Nsub, maskGlobal, pTxtGlobal, ...
        statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac);

    if saveFigures
        outPng = fullfile(grandOutFolder, sprintf('STATS_%s_meanSEM_byType_sig.png', thisExp));
        exportgraphics(fig2, outPng, 'BackgroundColor','white', 'Resolution',300);
        close(fig2);
        fprintf('[%s] Saved: %s\n', thisExp, outPng);
    end

    % --- Figure 3: Local+Global (2-panel) + both sig bars ---
fig3 = plot_local_global_onefig_with_sig( ...
    relTimes_ref, mu_sound, sem_sound, Ntr_local, ...
    mu_eff, sem_eff, Ntr_global, ...
    thisExp, Nsub, maskLocal, maskGlobal, maskInt, pTxtLocal, pTxtGlobal, ...
    statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac, thisExp);

    if saveFigures
        outPng = fullfile(grandOutFolder, sprintf('STATS_%s_meanSEM_LocalGlobal_sig.png', thisExp));
        lg_print_png_fixed(fig3, outPng, 300);
        close(fig3);
        fprintf('[%s] Saved: %s\n', thisExp, outPng);
    end



    % ===================== (A) CALLS TO ADD IN YOUR MAIN LOOP =====================
    % Put these calls INSIDE your for eidx=1:2 loop, AFTER you create fig2/fig3,
    % so we can store the needed data for E1 and E2 and then plot after the loop.

    % ----- ADD these lines right after you compute pTxtLocal / pTxtGlobal -----
    % (i.e., before plotting)
    if ~exist('COMPARE','var') || isempty(COMPARE)
        COMPARE = struct();
    end

    % Store E1/E2 TYPE data for later side-by-side comparison
    COMPARE.(thisExp).typeOrder     = typeOrder;
    COMPARE.(thisExp).mu_type       = mu_type;
    COMPARE.(thisExp).sem_type      = sem_type;
    COMPARE.(thisExp).maskGlobal    = maskGlobal;
    COMPARE.(thisExp).pTxtGlobal    = pTxtGlobal;
    COMPARE.(thisExp).maskInteract  = maskInt;    % --> Interaction effect
    COMPARE.(thisExp).Nsub          = Nsub;

    % Store E1/E2 GLOBAL-bySound-style (GS vs GD) data for standalone plot
    COMPARE.(thisExp).mu_eff        = mu_eff;
    COMPARE.(thisExp).sem_eff       = sem_eff;
    COMPARE.(thisExp).Ntr_global    = Ntr_global;


    % ===================== STORE FOR COMPARE FIGURE (UPDATED) =====================
    if ~exist('COMPARE','var') || isempty(COMPARE)
        COMPARE = struct();
    end

    % Store TYPE curves (the 4 LG conditions)
    COMPARE.(thisExp).typeOrder   = typeOrder;
    COMPARE.(thisExp).mu_type     = mu_type;
    COMPARE.(thisExp).sem_type    = sem_type;

    % Store LOCAL + GLOBAL significance (for top/bottom bars)
    COMPARE.(thisExp).maskLocal   = maskLocal;
    COMPARE.(thisExp).pTxtLocal   = pTxtLocal;

    COMPARE.(thisExp).maskGlobal  = maskGlobal;
    COMPARE.(thisExp).pTxtGlobal  = pTxtGlobal;

    % Meta
    COMPARE.(thisExp).Nsub        = Nsub;

    % (keep these if you still want your global-only standalone figure)
    COMPARE.(thisExp).mu_eff      = mu_eff;
    COMPARE.(thisExp).sem_eff     = sem_eff;
    COMPARE.(thisExp).Ntr_global  = Ntr_global;


end


% ----- ADD this call AFTER the for eidx=1:2 loop finishes -----
% (i.e., after the loop, before "DONE" print or after it)

figTypeCompare = plot_type_E1_vs_E2_side_by_side( ...
    COMPARE, relTimes_ref, statsWinMS, zoomPlotsToStatsWin, ...
    headroomFrac, barFracFromTop, pTxtYOffsetFrac, pText);

if saveFigures
    outPng = fullfile(grandOutFolder, 'STATS_E1E2_meanSEM_byType_sig_COMPARE.png');
    lg_print_png_fixed(figTypeCompare, outPng, 300);
    close(figTypeCompare);
    fprintf('[COMPARE] Saved: %s\n', outPng);
end

% Standalone global-only (GS vs GD) for E1
figE1GlobalOnly = plot_grand_mean_sem_global_only_with_sig( ...
    relTimes_ref, COMPARE.E1.mu_eff, COMPARE.E1.sem_eff, ...
    'E1', COMPARE.E1.Nsub, COMPARE.E1.Ntr_global, ...
    COMPARE.E1.maskGlobal, COMPARE.E1.maskInteract, COMPARE.E1.pTxtGlobal, ...
    statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac);

if saveFigures
    outPng = fullfile(grandOutFolder, 'STATS_E1_meanSEM_byGlobal_sig.png');
    exportgraphics(figE1GlobalOnly, outPng, 'BackgroundColor','white', 'Resolution',300);
    close(figE1GlobalOnly);
    fprintf('[E1] Saved: %s\n', outPng);
end

% Standalone global-only (GS vs GD) for E2 (optional but usually useful)
figE2GlobalOnly = plot_grand_mean_sem_global_only_with_sig( ...
    relTimes_ref, COMPARE.E2.mu_eff, COMPARE.E2.sem_eff, ...
    'E2', COMPARE.E2.Nsub, COMPARE.E2.Ntr_global, ...
    COMPARE.E2.maskGlobal, COMPARE.E2.maskInteract, COMPARE.E2.pTxtGlobal, ...
    statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac);

if saveFigures
    outPng = fullfile(grandOutFolder, 'STATS_E2_meanSEM_byGlobal_sig.png');
    exportgraphics(figE2GlobalOnly, outPng, 'BackgroundColor','white', 'Resolution',300);
    close(figE2GlobalOnly);
    fprintf('[E2] Saved: %s\n', outPng);
end


% --- 5 extra publication-style "zoomed-out" figures (NO CROPPING) ---%
%   - PUPIL_Bekinschtein_SSI_MATRIX.png
%   - PUPIL_BinnedPvals_SubjectMatrix.png
%   - PUPIL_Bekinschtein_SSI_MEAN_INDEX.png
pupil_export_5_extra_zoomed_figures( ...
    outpath, grandOutFolder, SSI, cfgSSI, cfgBIN, excludeSubs, string(useField));


% ===== OUTPUT FIGS (Bekinschtein-like) =====
cfgSSI.displayHideBinsAboveP = inf; % set e.g. 0.01 to hide weaker-than-0.01 bins

outMatPng = fullfile(grandOutFolder, 'PUPIL_Bekinschtein_SSI_MATRIX.png');
plot_pupil_subjectwise_matrix_post0(SSI, cfgSSI, outMatPng);

outIdxPng = fullfile(grandOutFolder, 'PUPIL_Bekinschtein_SSI_MEAN_INDEX.png');
plot_pupil_mean_subjectwise_index_post0(SSI, cfgSSI, outIdxPng);



% Binned matrix
if doBinnedMatrix
    pupil_export_binned_subject_matrix(string(outpath), string(grandOutFolder), cfgBIN, excludeSubs, string(useField));
end

% ---------- SAVE FULL TIME-BY-TIME INTERACTION t/p VECTORS ----------
interaction_timewise = struct();

expList = fieldnames(STATS);
for i = 1:numel(expList)
    expKey = expList{i};

    if isfield(STATS.(expKey), 'interactionSummary')
        interaction_timewise.(expKey).relTimes      = STATS.(expKey).relTimes;
        interaction_timewise.(expKey).statsWinMS    = STATS.(expKey).statsWinMS;
        interaction_timewise.(expKey).t_interaction = STATS.(expKey).interactionSummary.t;
        interaction_timewise.(expKey).p_interaction = STATS.(expKey).interactionSummary.p;
        interaction_timewise.(expKey).maskInteract  = STATS.(expKey).maskInteract;
        interaction_timewise.(expKey).clustersInteract = STATS.(expKey).clustersInteract;
        interaction_timewise.(expKey).sentence      = STATS.(expKey).interactionSummary.sentence;
    end
end

save(fullfile(grandOutFolder, 'STATS_PUPIL_interaction_timewise.mat'), ...
    'interaction_timewise', '-v7.3');


figInteraction = figure;
plot(interaction_timewise.E1.relTimes, interaction_timewise.E1.t_interaction, 'LineWidth', 2, 'Color',[0 0.45 0.85]);
hold on
plot(interaction_timewise.E2.relTimes, interaction_timewise.E2.t_interaction, 'LineWidth', 2, 'Color',[0.90 0.45 0.10]);
xline(0,'k-');
xlabel('Time (ms)');
ylabel('t value');
title('Interaction t-values across time');
grid on;
legend ('Counting', 'Visual Interference');

    if saveFigures
        outPngInteraction = fullfile(grandOutFolder, 'STATS_Interaction.png');
        exportgraphics(figInteraction, outPngInteraction, 'BackgroundColor','white', 'Resolution',300);
        close(figInteraction);
        fprintf('Figure Interaction Saved!');
    end


% Dilation comparison (Counting vs Visual Interference)
OUT_PUPIL_RAINBOX = pupil_export_rainbox_compare_E1E2( ...
    outpath, grandOutFolder, STATS, excludeSubs, useField);

% Separate 2 x 2 repeated-measures ANOVA
OUT_PUPIL_RAINBOX_ANOVA = pupil_export_rainbox_compare_E1E2_ANOVA( ...
    OUT_PUPIL_RAINBOX, grandOutFolder);


%% Summary rainbox

fprintf('\n=== CURRENT PUPIL RAINBOX DESCRIPTIVE SUMMARY ===\n');

T = OUT_PUPIL_RAINBOX.subjectTable;

pairs = {
    'Local Standard latency',      T.LocalStandard_Latency_E1_ms,  T.LocalStandard_Latency_E2_ms;
    'Local Deviant latency',       T.LocalDeviant_Latency_E1_ms,   T.LocalDeviant_Latency_E2_ms;
    'Local contrast latency',      T.LocalDeviant_Latency_E1_ms - T.LocalStandard_Latency_E1_ms, ...
                                   T.LocalDeviant_Latency_E2_ms - T.LocalStandard_Latency_E2_ms;

    'Local Standard mean',         T.LocalStandard_Mean_E1,        T.LocalStandard_Mean_E2;
    'Local Deviant mean',          T.LocalDeviant_Mean_E1,         T.LocalDeviant_Mean_E2;
    'Local contrast mean',         T.LocalDeviant_Mean_E1 - T.LocalStandard_Mean_E1, ...
                                   T.LocalDeviant_Mean_E2 - T.LocalStandard_Mean_E2;

    'Global Standard latency',     T.GlobalStandard_Latency_E1_ms, T.GlobalStandard_Latency_E2_ms;
    'Global Deviant latency',      T.GlobalDeviant_Latency_E1_ms,  T.GlobalDeviant_Latency_E2_ms;
    'Global contrast latency',     T.GlobalDeviant_Latency_E1_ms - T.GlobalStandard_Latency_E1_ms, ...
                                   T.GlobalDeviant_Latency_E2_ms - T.GlobalStandard_Latency_E2_ms;

    'Global Standard mean',        T.GlobalStandard_Mean_E1,       T.GlobalStandard_Mean_E2;
    'Global Deviant mean',         T.GlobalDeviant_Mean_E1,        T.GlobalDeviant_Mean_E2;
    'Global contrast mean',        T.GlobalDeviant_Mean_E1 - T.GlobalStandard_Mean_E1, ...
                                   T.GlobalDeviant_Mean_E2 - T.GlobalStandard_Mean_E2;
};

for i = 1:size(pairs,1)
    label = pairs{i,1};
    x1 = pairs{i,2};
    x2 = pairs{i,3};

    m = isfinite(x1) & isfinite(x2);
    x1 = x1(m);
    x2 = x2(m);

    n = numel(x1);
    if n < 2
        fprintf('%s | insufficient data\n', label);
        continue;
    end

    mu1 = mean(x1);
    mu2 = mean(x2);
    se1 = std(x1,0) / sqrt(n);
    se2 = std(x2,0) / sqrt(n);

    d = x1 - x2;
    dz = mean(d) / std(d,0);

    fprintf('%s | E1 = %.4f ± %.4f | E2 = %.4f ± %.4f | n = %d | dz = %.2f\n', ...
        label, mu1, se1, mu2, se2, n, dz);
end



%% Summaries +

% Statistics done
fprintf('\nDONE. Stats + outputs in: %s\n', grandOutFolder);
fprintf('\n=== INTERACTION CHECK ===\n');
fprintf('E1: significant samples = %d\n', nnz(STATS.E1.maskInteract));
fprintf('E2: significant samples = %d\n', nnz(STATS.E2.maskInteract));


%% SUMMARY RESULTS

% ===================== PUPIL RESULTS SUMMARY FOR PAPER =====================
fprintf('\n\n=== PUPIL MAIN SUMMARY ===\n');

for expName = ["E1","E2"]

    if expName == "E1"
        expLabel = "Counting";
    else
        expLabel = "Visual interference";
    end

    fprintf('\n%s\n', expLabel);
    fprintf('N subjects: %d\n', GRAND.(expName).Nsub);
    fprintf('Local trials: LS=%d | LD=%d\n', GRAND.(expName).Ntr_local(1), GRAND.(expName).Ntr_local(2));
    fprintf('Global trials: GS=%d | GD=%d\n', GRAND.(expName).Ntr_global(1), GRAND.(expName).Ntr_global(2));

    t = STATS.(expName).relTimes(:)';

    % ---------- LOCAL ----------
    maskL = STATS.(expName).maskLocal(:)';
    clsL  = STATS.(expName).clustersLocal;
    dLoc  = GRAND.(expName).mu_sound(2,:) - GRAND.(expName).mu_sound(1,:); % LD - LS

    fprintf('\n  Local effect:\n');
    if any(maskL)
        d = diff([false maskL false]);
        rs = find(d==1);
        re = find(d==-1)-1;

        for k = 1:numel(rs)
            fprintf('    Window %d: %.0f to %.0f ms\n', k, t(rs(k)), t(re(k)));
        end

        pAll = [clsL.pFWE];
        [pBest, iBest] = min(pAll);
        fprintf('    Best cluster pFWE = %.5f\n', pBest);

        idx = find(maskL);
        [~,ix] = max(abs(dLoc(idx)));
        idxPeak = idx(ix);
        fprintf('    Peak |LD-LS| in significant window = %.4f at %.0f ms\n', dLoc(idxPeak), t(idxPeak));

        if dLoc(idxPeak) > 0
            fprintf('    Direction: Local Deviant > Local Standard\n');
        else
            fprintf('    Direction: Local Standard > Local Deviant\n');
        end
    else
        fprintf('    No significant local cluster\n');
    end

    % ---------- GLOBAL ----------
    maskG = STATS.(expName).maskGlobal(:)';
    clsG  = STATS.(expName).clustersGlobal;
    dGlo  = GRAND.(expName).mu_effect(2,:) - GRAND.(expName).mu_effect(1,:); % GD - GS

    fprintf('\n  Global effect:\n');
    if any(maskG)
        d = diff([false maskG false]);
        rs = find(d==1);
        re = find(d==-1)-1;

        for k = 1:numel(rs)
            fprintf('    Window %d: %.0f to %.0f ms\n', k, t(rs(k)), t(re(k)));
        end

        pAll = [clsG.pFWE];
        [pBest, iBest] = min(pAll);
        fprintf('    Best cluster pFWE = %.5f\n', pBest);

        idx = find(maskG);
        [~,ix] = max(abs(dGlo(idx)));
        idxPeak = idx(ix);
        fprintf('    Peak |GD-GS| in significant window = %.4f at %.0f ms\n', dGlo(idxPeak), t(idxPeak));

        if dGlo(idxPeak) > 0
            fprintf('    Direction: Global Deviant > Global Standard\n');
        else
            fprintf('    Direction: Global Standard > Global Deviant\n');
        end
    else
        fprintf('    No significant global cluster\n');
    end

    % ---------- INTERACTION ----------
    fprintf('\n  Interaction:\n');
    if isfield(STATS.(expName), 'interactionSummary')
        disp(STATS.(expName).interactionSummary.table)
    else
        fprintf('    No interaction summary found\n');
    end
end


% ===================== PUPIL INDIVIDUAL-LEVEL SUMMARY =====================

fprintf('\n\n=== PUPIL INDIVIDUAL-LEVEL SUMMARY ===\n');

for expName = ["E1","E2"]

    if expName == "E1"
        expLabel = "Counting";
    else
        expLabel = "Visual interference";
    end

    t = SSI.(expName).times(:)';
    if all(isnan(t))
        continue;
    end

    % usually post-0 only for SSI interpretation
    keep = t >= 0;
    t2 = t(keep);

    dt = median(diff(t2)); % ms per sample

    L = SSI.(expName).local(:,keep);
    G = SSI.(expName).global(:,keep);

    hasLoc = any(L > 0, 2);
    hasGlo = any(G > 0, 2);

    durLoc = sum(L > 0, 2) * dt;
    durGlo = sum(G > 0, 2) * dt;

    meanLoc = mean(L, 1, 'omitnan');
    meanGlo = mean(G, 1, 'omitnan');

    [peakLoc, ixLoc] = max(meanLoc);
    [peakGlo, ixGlo] = max(meanGlo);

    fprintf('\n%s\n', expLabel);
    fprintf('  Local | subjects with any effect: %d/%d (%.1f%%)\n', ...
        sum(hasLoc), size(L,1), 100*sum(hasLoc)/size(L,1));
    fprintf('  Local | peak mean SSI: %.3f at %.0f ms\n', peakLoc, t2(ixLoc));
    fprintf('  Local | median duration among all subjects: %.0f ms\n', median(durLoc,'omitnan'));
    fprintf('  Local | median duration among subjects with effect: %.0f ms\n', median(durLoc(hasLoc),'omitnan'));

    fprintf('  Global | subjects with any effect: %d/%d (%.1f%%)\n', ...
        sum(hasGlo), size(G,1), 100*sum(hasGlo)/size(G,1));
    fprintf('  Global | peak mean SSI: %.3f at %.0f ms\n', peakGlo, t2(ixGlo));
    fprintf('  Global | median duration among all subjects: %.0f ms\n', median(durGlo,'omitnan'));
    fprintf('  Global | median duration among subjects with effect: %.0f ms\n', median(durGlo(hasGlo),'omitnan'));
end




%% Grand average summary

% Dilation comparison (Counting vs Visual Interference)
OUT_PUPIL_RAINBOX = pupil_export_rainbox_compare_E1E2( ...
    outpath, grandOutFolder, STATS, excludeSubs, useField);

fprintf('\n=== CURRENT PUPIL RAINBOX DESCRIPTIVE SUMMARY ===\n');

T = OUT_PUPIL_RAINBOX.subjectTable;

pairs = {
    'Local Standard latency',      T.LocalStandard_Latency_E1_ms,  T.LocalStandard_Latency_E2_ms;
    'Local Deviant latency',       T.LocalDeviant_Latency_E1_ms,   T.LocalDeviant_Latency_E2_ms;
    'Local contrast latency',      T.LocalDeviant_Latency_E1_ms - T.LocalStandard_Latency_E1_ms, ...
                                   T.LocalDeviant_Latency_E2_ms - T.LocalStandard_Latency_E2_ms;

    'Local Standard mean',         T.LocalStandard_Mean_E1,        T.LocalStandard_Mean_E2;
    'Local Deviant mean',          T.LocalDeviant_Mean_E1,         T.LocalDeviant_Mean_E2;
    'Local contrast mean',         T.LocalDeviant_Mean_E1 - T.LocalStandard_Mean_E1, ...
                                   T.LocalDeviant_Mean_E2 - T.LocalStandard_Mean_E2;

    'Global Standard latency',     T.GlobalStandard_Latency_E1_ms, T.GlobalStandard_Latency_E2_ms;
    'Global Deviant latency',      T.GlobalDeviant_Latency_E1_ms,  T.GlobalDeviant_Latency_E2_ms;
    'Global contrast latency',     T.GlobalDeviant_Latency_E1_ms - T.GlobalStandard_Latency_E1_ms, ...
                                   T.GlobalDeviant_Latency_E2_ms - T.GlobalStandard_Latency_E2_ms;

    'Global Standard mean',        T.GlobalStandard_Mean_E1,       T.GlobalStandard_Mean_E2;
    'Global Deviant mean',         T.GlobalDeviant_Mean_E1,        T.GlobalDeviant_Mean_E2;
    'Global contrast mean',        T.GlobalDeviant_Mean_E1 - T.GlobalStandard_Mean_E1, ...
                                   T.GlobalDeviant_Mean_E2 - T.GlobalStandard_Mean_E2;
};

for i = 1:size(pairs,1)
    label = pairs{i,1};
    x1 = pairs{i,2};
    x2 = pairs{i,3};

    m = isfinite(x1) & isfinite(x2);
    x1 = x1(m);
    x2 = x2(m);

    n = numel(x1);
    if n < 2
        fprintf('%s | insufficient data\n', label);
        continue;
    end

    mu1 = mean(x1);
    mu2 = mean(x2);
    se1 = std(x1,0) / sqrt(n);
    se2 = std(x2,0) / sqrt(n);

    d = x1 - x2;
    dz = mean(d) / std(d,0);

    fprintf('%s | E1 = %.4f ± %.4f | E2 = %.4f ± %.4f | n = %d | dz = %.2f\n', ...
        label, mu1, se1, mu2, se2, n, dz);
end

fprintf('\nDONE. Stats + outputs in: %s\n', grandOutFolder);
fprintf('\n=== INTERACTION CHECK ===\n');
fprintf('E1: significant samples = %d\n', nnz(STATS.E1.maskInteract));
fprintf('E2: significant samples = %d\n', nnz(STATS.E2.maskInteract));

% ---------- EXPORT ONE SINGLE INTERACTION TABLE ----------
Tinteraction = export_interaction_summary_table(STATS, grandOutFolder);
disp(Tinteraction);

% NEW: richer descriptive summaries of significant group clusters
PUPIL_GROUP_DESC = pupil_export_group_cluster_descriptives( ...
    outpath, grandOutFolder, STATS, excludeSubs, useField);

% NEW: supplementary tables
SUPP_PUPIL = pupil_export_supplementary_tables_all( ...
    outpath, grandOutFolder, STATS, GRAND, SSI, Tinteraction, ...
    OUT_PUPIL_RAINBOX, PUPIL_GROUP_DESC, excludeSubs);


%% ===================== FUNCTIONS =====================

function m = subject_mean_curve(T, fieldName, condVar, condValue)
% subject_mean_curve
%   Computes a subject-level mean curve (1 x nTime) for trials matching condValue.
%   If the subject has no trials for that condition -> returns NaNs.

condVec = string(T.(condVar));
idx = (condVec == string(condValue));

nTime = numel(T.(fieldName){1});
if ~any(idx)
    m = nan(1, nTime);
    return;
end

trialCells = T.(fieldName)(idx);
pupMat = vertcat(trialCells{:}); % [nTrials x nTime]
m = mean(pupMat, 1, 'omitnan');
end

function [mu, sem] = mean_sem_over_subjects(matsCell)
% mean_sem_over_subjects
% matsCell: {M1, M2, ...} each Mi is [nSub x nTime]
% Returns:
%   mu  = [nCond x nTime] mean across subjects
%   sem = [nCond x nTime] SEM across subjects (time-resolved)

nCond = numel(matsCell);
nTime = size(matsCell{1},2);

mu  = nan(nCond, nTime);
sem = nan(nCond, nTime);

for c = 1:nCond
    M = matsCell{c};
    mu(c,:) = mean(M, 1, 'omitnan');

    n_t  = sum(~isnan(M), 1);
    sd_t = std(M, 0, 1, 'omitnan');
    sem(c,:) = sd_t ./ sqrt(max(n_t,1));
end
end

%% ---------- CLUSTER PERMUTATION TESTS (1D across time) ----------
function [sigMask, clusters] = clusterperm_paired_1d(A, B, alphaCluster, alphaFWE, nPerm)
% paired test implemented as one-sample on D=A-B
D = A - B;
[sigMask, clusters] = clusterperm_onesample_1d(D, alphaCluster, alphaFWE, nPerm);
end

function [sigMask, clusters] = clusterperm_onesample_1d(D, alphaCluster, alphaFWE, nPerm)
% clusterperm_onesample_1d
%   One-sample cluster permutation test against 0, across time (columns).
%   Uses sign-flip permutations, cluster mass = sum(|t|) in cluster.

allnan = all(isnan(D),2);
D = D(~allnan,:);
nSub  = size(D,1);
nTime = size(D,2);

sigMask = false(1,nTime);
clusters = struct('startIdx',{},'endIdx',{},'clusterMass',{},'pFWE',{});

if nSub < 2
    warning('clusterperm_onesample_1d: Not enough subjects (n=%d). Returning empty mask.', nSub);
    return;
end

tObs = tstat_onesample(D);

df = nSub - 1;
tcrit = tinv(1 - alphaCluster/2, df);

supra = abs(tObs) > tcrit;
obsClusters = find_clusters_1d(supra);
if isempty(obsClusters), return; end

obsMass = zeros(numel(obsClusters),1);
for i = 1:numel(obsClusters)
    idx = obsClusters(i).startIdx:obsClusters(i).endIdx;
    obsMass(i) = sum(abs(tObs(idx)));
end

maxMass = zeros(nPerm,1);
for p = 1:nPerm
    flips = (rand(nSub,1) > 0.5);
    sgn = ones(nSub,1);
    sgn(flips) = -1;

    Dp = D .* sgn;
    tP = tstat_onesample(Dp);

    supraP = abs(tP) > tcrit;
    clP = find_clusters_1d(supraP);

    if isempty(clP)
        maxMass(p) = 0;
    else
        masses = zeros(numel(clP),1);
        for k = 1:numel(clP)
            idx = clP(k).startIdx:clP(k).endIdx;
            masses(k) = sum(abs(tP(idx)));
        end
        maxMass(p) = max(masses);
    end
end

pFWE = zeros(numel(obsClusters),1);
for i = 1:numel(obsClusters)
    pFWE(i) = mean(maxMass >= obsMass(i));
end

keep = pFWE < alphaFWE;

for i = 1:numel(obsClusters)
    clusters(i).startIdx    = obsClusters(i).startIdx;
    clusters(i).endIdx      = obsClusters(i).endIdx;
    clusters(i).clusterMass = obsMass(i);
    clusters(i).pFWE        = pFWE(i);
end

for i = find(keep(:))'
    sigMask(obsClusters(i).startIdx:obsClusters(i).endIdx) = true;
end
end

function t = tstat_onesample(D)
% time-wise one-sample t-stat against 0
n  = sum(~isnan(D),1);
mu = mean(D,1,'omitnan');
sd = std(D,0,1,'omitnan');

t = mu ./ (sd ./ sqrt(max(n,1)));
t(n < 2) = 0;
end

function clusters = find_clusters_1d(mask)
% find contiguous true segments
mask = mask(:)';
d = diff([false mask false]);
starts = find(d==1);
ends   = find(d==-1)-1;
clusters = struct('startIdx',num2cell(starts), 'endIdx',num2cell(ends));
end

%% ---------- WINDOWED WRAPPERS (NEW) ----------
function [sigMaskFull, clustersFull] = clusterperm_paired_1d_windowed(A, B, timeVector, statsWinMS, alphaCluster, alphaFWE, nPerm)
% Runs paired clusterperm only within statsWinMS, maps mask back to full axis.

idxWin = (timeVector >= statsWinMS(1)) & (timeVector <= statsWinMS(2));
if ~any(idxWin)
    error('clusterperm_paired_1d_windowed: statsWinMS [%g %g] has no samples in timeVector.', statsWinMS(1), statsWinMS(2));
end

Awin = A(:, idxWin);
Bwin = B(:, idxWin);

[sigMaskWin, clustersWin] = clusterperm_paired_1d(Awin, Bwin, alphaCluster, alphaFWE, nPerm);

sigMaskFull = false(1, numel(timeVector));
sigMaskFull(idxWin) = sigMaskWin;

clustersFull = remap_clusters_to_full(clustersWin, idxWin);
end

function [sigMaskFull, clustersFull] = clusterperm_onesample_1d_windowed(D, timeVector, statsWinMS, alphaCluster, alphaFWE, nPerm)
% Runs one-sample clusterperm only within statsWinMS, maps mask back to full axis.

idxWin = (timeVector >= statsWinMS(1)) & (timeVector <= statsWinMS(2));
if ~any(idxWin)
    error('clusterperm_onesample_1d_windowed: statsWinMS [%g %g] has no samples in timeVector.', statsWinMS(1), statsWinMS(2));
end

Dwin = D(:, idxWin);

[sigMaskWin, clustersWin] = clusterperm_onesample_1d(Dwin, alphaCluster, alphaFWE, nPerm);

sigMaskFull = false(1, numel(timeVector));
sigMaskFull(idxWin) = sigMaskWin;

clustersFull = remap_clusters_to_full(clustersWin, idxWin);
end

function clustersFull = remap_clusters_to_full(clustersWin, idxWin)
% Converts cluster start/end indices from WINDOW indexing to FULL indexing.

clustersFull = clustersWin;
if isempty(clustersWin)
    return;
end

fullIdx = find(idxWin); % mapping: winIndex -> fullIndex
for i = 1:numel(clustersWin)
    clustersFull(i).startIdx = fullIdx(clustersWin(i).startIdx);
    clustersFull(i).endIdx   = fullIdx(clustersWin(i).endIdx);
end
end

%% ---------- P-TEXT (ONLY "pFWE < 0.05") ----------
function pTxt = ptext_simple(clusters, alphaFWE)
pTxt = "";
if isempty(clusters), return; end
pVals = [clusters.pFWE];
pVals = pVals(~isnan(pVals));
if isempty(pVals), return; end
if min(pVals) < alphaFWE
    pTxt = sprintf('p_{FWE} < %.2f', alphaFWE);
end
end

%% ---------- PLOTTING ----------
function fig = plot_grand_mean_sem_local_with_sig(timeVector, soundOrder, grandMean, grandSEM, expLabel, Nsub, Ntr_local, sigMaskLocal, sigMaskInteract, pTxtLocal, statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac)
% Local plot: LS vs LD with local significance bar

% colors = [0 0.4470 0.7410;
%     0.8500 0.3250 0.0980];

colors = [0 0.45 0.85;
    0.90 0.45 0.10];

fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w', 'InvertHardcopy','off');
ax = axes(fig); hold(ax,'on');
set(ax,'Color','w','XColor','k','YColor','k');
set(groot,'defaultAxesFontSize',11);

title(ax, sprintf("Grand average | %s | Local | Nsub=%d | Ntrials LS=%d LD=%d", ...
    expLabel, Nsub, Ntr_local(1), Ntr_local(2)), FontSize=12);

xlabel(ax,"Time (ms)", FontSize=12);
ylabel(ax,"Pupil size (z-score)", FontSize=12);

lineH = gobjects(numel(soundOrder),1);

for s = 1:numel(soundOrder)
    mu  = grandMean(s, :);
    sem = grandSEM(s, :);
    if all(isnan(mu)), continue; end

    upper = mu + sem;
    lower = mu - sem;

    hFill = fill(ax, [timeVector fliplr(timeVector)], [upper fliplr(lower)], ...
        colors(s,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    set(hFill,'HandleVisibility','off');

    lineH(s) = plot(ax, timeVector, mu, 'LineWidth', 2, 'Color', colors(s,:), ...
        'DisplayName', char(soundOrder(s)));
end

xline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');

% add the top bar (this function also adds ylim headroom automatically)
add_sig_bar_top(ax, timeVector, sigMaskLocal, headroomFrac, barFracFromTop, [0 0 1], pTxtLocal, pTxtYOffsetFrac);
% Interaction bar (gray) shown BELOW the blue bar if significant
add_sig_bar_top_secondary(ax, timeVector, sigMaskInteract, barFracFromTop, 0.045, [0.50 0.50 0.50]);

% --- FIXED Y-SCALE FOR MANUSCRIPT FIGURES ---
ylim(ax, [-0.30 0.70]);

if zoomPlotsToStatsWin
    xlim(ax, statsWinMS);
end

legend(ax, lineH(isgraphics(lineH)), 'Location','southwest', FontSize=12);
grid(ax,'on');
box(ax,'off');              % removes top/right border lines
ax.TickDir = 'out';         % (optional) matches a cleaner style
ax.Layer   = 'top';         % keeps grid behind data nicely

hold(ax,'off');
end

function fig = plot_grand_mean_sem_by_type_with_sig(timeVector, typeOrder, grandMean, grandSEM, expLabel, Nsub, sigMaskGlobal, pTxtGlobal, statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac)
% Type plot: 4 lines (FLS/RLS/FLD/RLD) with a GLOBAL (GS vs GD) significance bar

% colors = [ ...
%     0.0000 0.4470 0.7410;
%     0.0000 0.4470 0.7410;
%     0.8500 0.3250 0.0980;
%     0.8500 0.3250 0.0980];

colors = [ ...
    0 0.45 0.85;
    0 0.45 0.85;
    0.90 0.45 0.10;
    0.90 0.45 0.10];

lineStyles = {'--','-','--','-'};

fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w', 'InvertHardcopy','off');
ax = axes(fig); hold(ax,'on');
set(ax,'Color','w','XColor','k','YColor','k');
set(groot,'defaultAxesFontSize',11);

title(ax, sprintf("Grand average | %s | Pupil response by type | Nsub=%d", expLabel, Nsub), FontSize=12);
xlabel(ax,"Time (ms)", FontSize=12);
ylabel(ax,"Pupil size (z-score)", FontSize=12);

lineH = gobjects(numel(typeOrder),1);

for k = 1:numel(typeOrder)
    mu  = grandMean(k, :);
    sem = grandSEM(k, :);
    if all(isnan(mu)), continue; end

    upper = mu + sem;
    lower = mu - sem;

    hFill = fill(ax, [timeVector fliplr(timeVector)], [upper fliplr(lower)], ...
        colors(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    set(hFill,'HandleVisibility','off');

    lineH(k) = plot(ax, timeVector, mu, ...
        'Color', colors(k,:), 'LineStyle', lineStyles{k}, 'LineWidth', 2, ...
        'DisplayName', char(typeOrder(k)));
end

xline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');

% IMPORTANT: the bar here indicates GLOBAL effect (GS vs GD), not "one line vs others"
add_sig_bar_top(ax, timeVector, sigMaskGlobal, headroomFrac, barFracFromTop, [0 0 1], pTxtGlobal, pTxtYOffsetFrac);

if zoomPlotsToStatsWin
    xlim(ax, statsWinMS);
end

legend(ax, lineH(isgraphics(lineH)), "Location","southwest", FontSize=12);
grid(ax,'on');
box(ax,'off');              % removes top/right border lines
ax.TickDir = 'out';         % (optional) matches a cleaner style
ax.Layer   = 'top';         % keeps grid behind data nicely

hold(ax,'off');
end


function [ssiLevel, pVec] = pupil_individual_pmap_bekinschtein(T, fieldName, timeVector, condVec, condA, condB, cfg)
% pupil_individual_pmap_bekinschtein
% Builds a Bekinschtein-like "significance index" time-series for ONE subject:
%   - sample-by-sample two-sample t-tests across TRIALS between condA and condB
%   - optional random trial balancing
%   - optional Student or Welch t-test
%   - extent rule: p<=alpha for >=minConsec consecutive samples
%   - binning: for each timepoint passing extent rule, assign the MOST stringent p-bin met
%
% cfg options:
%   cfg.testMethod    = "welch" or "student"
%   cfg.balanceTrials = true or false
%   cfg.balanceMode   = "elementary" or "contrast"
%
% balanceMode:
%   "contrast"   -> equalizes only condA vs condB trial numbers.
%   "elementary" -> equalizes FLS, RLS, FLD, RLD first, then builds the contrast.
%                   This is closer to the factorial LG/Bekinschtein-style logic.

arguments
    T table
    fieldName (1,1) string
    timeVector (1,:) double
    condVec
    condA (1,1) string
    condB (1,1) string
    cfg struct
end

% ---------- defaults ----------
if ~isfield(cfg,'alpha_ind'),      cfg.alpha_ind = 0.05; end
if ~isfield(cfg,'minConsec_ind'),  cfg.minConsec_ind = 10; end
if ~isfield(cfg,'pBins_matrix'),   cfg.pBins_matrix = [0.05 0.01 0.005 0.001 0.0005 0.0001]; end
if ~isfield(cfg,'post0_only'),     cfg.post0_only = true; end

if ~isfield(cfg,'testMethod') || isempty(cfg.testMethod)
    cfg.testMethod = "welch";
end
if ~isfield(cfg,'balanceTrials') || isempty(cfg.balanceTrials)
    cfg.balanceTrials = false;
end
if ~isfield(cfg,'balanceMode') || isempty(cfg.balanceMode)
    cfg.balanceMode = "contrast";
end

condVec = string(condVec);

nTime = numel(timeVector);
pVec  = nan(1,nTime);
ssiLevel = zeros(1,nTime);

% ---------- select trials, optionally balanced ----------
[idxA, idxB] = pupil_get_balanced_trial_indices(T, condVec, condA, condB, cfg);

if ~any(idxA) || ~any(idxB)
    return;
end

XA = vertcat(T.(fieldName){idxA}); % [nA x nTime]
XB = vertcat(T.(fieldName){idxB}); % [nB x nTime]

% Safety: ensure correct nTime
if size(XA,2) ~= nTime || size(XB,2) ~= nTime
    error('Trial vectors in %s do not match timeVector length.', fieldName);
end

% ---------- choose Student vs Welch ----------
switch lower(string(cfg.testMethod))
    case {"student","pooled","equal"}
        vartype = 'equal';    % classic Student two-sample t-test
    case {"welch","unequal"}
        vartype = 'unequal';  % Welch two-sample t-test
    otherwise
        error('Unknown cfg.testMethod = %s. Use "welch" or "student".', string(cfg.testMethod));
end

% ---------- sample-by-sample two-sample t-test across trials ----------
for tt = 1:nTime
    a = XA(:,tt);
    b = XB(:,tt);

    a = a(~isnan(a));
    b = b(~isnan(b));

    if numel(a) < 2 || numel(b) < 2
        pVec(tt) = NaN;
    else
        [~,pVec(tt)] = ttest2(a, b, 'Vartype', vartype);
    end
end

% ---------- apply extent rule on FIRST threshold ----------
alpha = cfg.alpha_ind;
mask = (pVec <= alpha);

% Optional: only enforce extent in post-0 region
if cfg.post0_only
    mask(timeVector < 0) = false;
end

maskExt = apply_extent_1d(mask, cfg.minConsec_ind);

% ---------- bin only where extent holds ----------
pb = cfg.pBins_matrix(:)'; 
nBins = numel(pb);

for tt = 1:nTime
    if ~maskExt(tt) || isnan(pVec(tt))
        continue;
    end

    lvl = 0;
    for k = 1:nBins
        if pVec(tt) <= pb(k)
            lvl = k;
        end
    end

    ssiLevel(tt) = lvl;
end

end


function [idxA, idxB] = pupil_get_balanced_trial_indices(T, condVec, condA, condB, cfg)
% Returns logical indices for condA and condB.
% If cfg.balanceTrials=true, trials are randomly subsampled.
%
% balanceMode = "contrast":
%   Equalizes only condA and condB trial counts.
%
% balanceMode = "elementary":
%   Equalizes FLS, RLS, FLD, RLD first.
%   This is useful for local/global contrasts because it preserves the
%   factorial local-global structure before testing.

condVec = string(condVec);
condA   = string(condA);
condB   = string(condB);

idxA0 = (condVec == condA);
idxB0 = (condVec == condB);

idxA = idxA0;
idxB = idxB0;

if ~isfield(cfg,'balanceTrials') || ~cfg.balanceTrials
    return;
end

mode = lower(string(cfg.balanceMode));

switch mode

    case "elementary"

        % This requires the four elementary LG labels in T.type
        if ~ismember('type', T.Properties.VariableNames)
            warning('balanceMode="elementary" requested, but T.type is missing. Falling back to contrast balancing.');
            [idxA, idxB] = pupil_balance_two_sets(idxA0, idxB0);
            return;
        end

        typeVec = string(T.type);

        elementaryTypes = [ ...
            "Frequent Local Standard", ...
            "Rare Local Standard", ...
            "Frequent Local Deviant", ...
            "Rare Local Deviant"];

        nEach = zeros(1,numel(elementaryTypes));
        for k = 1:numel(elementaryTypes)
            nEach(k) = sum(typeVec == elementaryTypes(k));
        end

        nMin = min(nEach);

        if nMin < 2
            warning('Not enough trials in one or more elementary LG cells. Falling back to contrast balancing.');
            [idxA, idxB] = pupil_balance_two_sets(idxA0, idxB0);
            return;
        end

        keepElementary = false(height(T),1);

        for k = 1:numel(elementaryTypes)
            idxK = find(typeVec == elementaryTypes(k));
            idxK = idxK(randperm(numel(idxK), nMin));
            keepElementary(idxK) = true;
        end

        idxA = idxA0 & keepElementary;
        idxB = idxB0 & keepElementary;

    case "contrast"

        [idxA, idxB] = pupil_balance_two_sets(idxA0, idxB0);

    otherwise

        error('Unknown cfg.balanceMode = %s. Use "elementary" or "contrast".', string(cfg.balanceMode));
end

end


function [idxA, idxB] = pupil_balance_two_sets(idxA0, idxB0)
% Randomly subsamples two condition sets to the same number of trials.

idxA = false(size(idxA0));
idxB = false(size(idxB0));

ia = find(idxA0);
ib = find(idxB0);

nMin = min(numel(ia), numel(ib));

if nMin < 2
    return;
end

iaKeep = ia(randperm(numel(ia), nMin));
ibKeep = ib(randperm(numel(ib), nMin));

idxA(iaKeep) = true;
idxB(ibKeep) = true;

end


function maskExt = apply_extent_1d(mask, minConsec)
% Keeps only segments of TRUE with length >= minConsec
mask = mask(:)';
d = diff([false mask false]);
starts = find(d==1);
ends   = find(d==-1)-1;

maskExt = false(size(mask));
for i = 1:numel(starts)
    if (ends(i)-starts(i)+1) >= minConsec
        maskExt(starts(i):ends(i)) = true;
    end
end
end

function cmap = pbin_cmap_simple(cfg, nColors)
% pbin_cmap_simple
% Uses user-defined bin colors if provided in cfg.binColors.
% nColors is expected to be nBins+1 (levels 0..nBins).

% Default fallback
if nargin < 2 || isempty(nColors)
    nColors = 7;
end
cmap = jet(max(nColors,2)); %#ok<JET>
cmap(1,:) = [0.10 0.45 0.90]; % default "none"

% Use custom colors if present and valid
if nargin >= 1 && isstruct(cfg) && isfield(cfg,'binColors') && ~isempty(cfg.binColors)
    C = cfg.binColors;
    if size(C,2) ~= 3
        error('cfg.binColors must be [nColors x 3] RGB.');
    end
    if size(C,1) ~= nColors
        error('cfg.binColors has %d rows but needs %d (nBins+1).', size(C,1), nColors);
    end
    cmap = C;
end
end

function add_two_row_colorbars_right_axes(axTopRight, axBotRight, pBins, cmap, labelStr, clim)
% Robust: create colorbars using AXES handles (not the FIG handle).
% Also includes level 0 ("none") in the tick labels.

nBins = numel(pBins);
if nargin < 6, clim = [-0.5 nBins+0.5]; end

% Tick positions represent levels: 0..nBins
ticks = 0:nBins;

% Labels: 0 = none, then <=pBins(k)
tickLabels = strings(size(ticks));
tickLabels(1) = ">0.05";
for k = 1:nBins
    tickLabels(k+1) = sprintf('<=%.4g', pBins(k));
end

% --- TOP colorbar (use top-right axes as anchor) ---
cb1 = colorbar(axTopRight);
cb1.Position = [0.92 0.56 0.015 0.35];
colormap(axTopRight, cmap);
cb1.Limits = clim;
cb1.Ticks  = ticks;
cb1.TickLabels = tickLabels;
ylabel(cb1, labelStr);

% --- BOTTOM colorbar (use bottom-right axes as anchor) ---
cb2 = colorbar(axBotRight);
cb2.Position = [0.92 0.10 0.015 0.35];
colormap(axBotRight, cmap);
cb2.Limits = clim;
cb2.Ticks  = ticks;
cb2.TickLabels = tickLabels;
ylabel(cb2, labelStr);
end

function plot_pupil_subjectwise_matrix_post0(SSI, cfg, outPng)
% 2x2 matrix: subjects x time, binned SSI levels

tms = SSI.E1.times;
keep = true(size(tms));
if isfield(cfg,'post0_only') && cfg.post0_only
    keep = tms >= 0;
end
tms2 = tms(keep);

E1L = SSI.E1.local(:,keep);
E1G = SSI.E1.global(:,keep);
E2L = SSI.E2.local(:,keep);
E2G = SSI.E2.global(:,keep);

% Optional: hide bins weaker than some threshold
if isfield(cfg,'displayHideBinsAboveP') && isfinite(cfg.displayHideBinsAboveP)
    pb = cfg.pBins_matrix(:)';
    hideLvl = find(pb <= cfg.displayHideBinsAboveP, 1, 'last');
    if isempty(hideLvl), hideLvl = 0; end
    E1L(E1L < hideLvl) = 0; E1G(E1G < hideLvl) = 0;
    E2L(E2L < hideLvl) = 0; E2G(E2G < hideLvl) = 0;
end

pb   = cfg.pBins_matrix(:)';   % e.g. [0.05 0.01 0.005 ...]
nBins = numel(pb);

% We want levels 0..nBins => need nBins+1 colors
cmap = pbin_cmap_simple(cfg, nBins+1);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
t.OuterPosition = [0.03 0.06 0.88 0.90]; % more top margin to avoid overlap

ax1 = nexttile;
imagesc(tms2, 1:size(E1L,1), E1L);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants');
title('Counting | Local effects', 'FontSize', 10, 'FontWeight','bold');
colormap(ax1,cmap); caxis(ax1,[-0.5 nBins+0.5]);

ax2 = nexttile;
imagesc(tms2, 1:size(E1G,1), E1G);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants');
title('Counting | Global effects', 'FontSize', 10, 'FontWeight','bold');
colormap(ax2,cmap); caxis(ax2,[-0.5 nBins+0.5]);

ax3 = nexttile;
imagesc(tms2, 1:size(E2L,1), E2L);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants');
title('Visual interference | Local effects', 'FontSize', 10, 'FontWeight','bold');
colormap(ax3,cmap); caxis(ax3,[-0.5 nBins+0.5]);

ax4 = nexttile;
imagesc(tms2, 1:size(E2G,1), E2G);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants');
title('Visual interference | Global effects', 'FontSize', 10, 'FontWeight','bold');
colormap(ax4,cmap); caxis(ax4,[-0.5 nBins+0.5]);

% Colorbars: attach to axes (robust across MATLAB versions)
lbl = 'p-values';
if isfield(cfg,'binLabel') && strlength(string(cfg.binLabel))>0
    lbl = char(cfg.binLabel);
end
add_two_row_colorbars_right_axes(ax2, ax4, pb, cmap, lbl, [-0.5 nBins+0.5]);

% sgtitle(sprintf('PUPIL Bekinschtein-like SUBJECT matrix (subjects×time) | minConsec=%d | post-0 only', ...
%     cfg.minConsec_ind), 'FontWeight','bold', 'FontSize', 12);

sgtitle('Local and global effects in individual subjects | Pupil size', 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function plot_pupil_mean_subjectwise_index_post0(SSI, cfg, outPng)
% Mean SSI index across subjects: Local and Global, for E1 and E2 panels

tms = SSI.E1.times;
keep = true(size(tms));
if isfield(cfg,'post0_only') && cfg.post0_only
    keep = tms >= 0;
end
tms2 = tms(keep);

mE1L = mean(SSI.E1.local(:,keep),1,'omitnan');
mE1G = mean(SSI.E1.global(:,keep),1,'omitnan');
mE2L = mean(SSI.E2.local(:,keep),1,'omitnan');
mE2G = mean(SSI.E2.global(:,keep),1,'omitnan');

yMax = max([mE1L(:); mE1G(:); mE2L(:); mE2G(:)]);
yMax = max(1, 1.05*yMax);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(tms2, mE1L, 'LineWidth', 2, 'Color',[0 0.45 0.85]);
plot(tms2, mE1G, 'LineWidth', 2, 'Color',[0.90 0.45 0.10]);
xlabel('Time after 5th sound onset (ms)');
ylabel('Mean binned index');
title('Counting');
legend({'Local effect','Global effect'}, 'Location','northeast');
ylim([0 yMax]);

nexttile; hold on; grid on;
plot(tms2, mE2L, 'LineWidth', 2, 'Color',[0 0.45 0.85]);
plot(tms2, mE2G, 'LineWidth', 2, 'Color',[0.90 0.45 0.10]);
xlabel('Time after 5th sound onset (ms)');
ylabel('Mean binned index');
title('Visual interference');
legend({'Local effect','Global effect'}, 'Location','northeast');
ylim([0 yMax]);

% sgtitle('PUPIL Mean subject-wise binned index | post-0 only', 'FontWeight','bold');
sgtitle('Mean individual statistical significance index | Pupil size', 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end


function fig = plot_grand_mean_sem_global_only_with_sig(timeVector, mu_eff, sem_eff, expLabel, Nsub, Ntr_global, sigMaskGlobal, sigMaskInteract, pTxtGlobal, statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac)
% plot_grand_mean_sem_global_only_with_sig
%   Standalone global-only plot (GS vs GD), same style as your Local-bySound figure,
%   but using the global effect curves (from Frequent vs Rare mapping).

% cStd = [0 0.4470 0.7410];
% cDev = [0.8500 0.3250 0.0980];
cStd = [0 0.45 0.85];
cDev = [0.90 0.45 0.10];


fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w', 'InvertHardcopy','off');
ax = axes(fig); hold(ax,'on');
set(ax,'Color','w','XColor','k','YColor','k');
set(groot,'defaultAxesFontSize',11);

title(ax, sprintf("Grand average | %s | Global | Nsub=%d | Ntrials GS=%d GD=%d", ...
    expLabel, Nsub, Ntr_global(1), Ntr_global(2)), FontSize=12);

xlabel(ax,"Time (ms)", FontSize=12);
ylabel(ax,"Pupil size (z-score)", FontSize=12);

% GS
upper = mu_eff(1,:) + sem_eff(1,:);
lower = mu_eff(1,:) - sem_eff(1,:);
hf = fill(ax,[timeVector fliplr(timeVector)],[upper fliplr(lower)],cStd,'FaceAlpha',0.15,'EdgeColor','none');
set(hf,'HandleVisibility','off');
hGS = plot(ax,timeVector,mu_eff(1,:),'Color',cStd,'LineWidth',2,'DisplayName','Global Standard');

% GD
upper = mu_eff(2,:) + sem_eff(2,:);
lower = mu_eff(2,:) - sem_eff(2,:);
hf = fill(ax,[timeVector fliplr(timeVector)],[upper fliplr(lower)],cDev,'FaceAlpha',0.15,'EdgeColor','none');
set(hf,'HandleVisibility','off');
hGD = plot(ax,timeVector,mu_eff(2,:),'Color',cDev,'LineWidth',2,'DisplayName','Global Deviant');

% Global bar = GS vs GD
add_sig_bar_top(ax, timeVector, sigMaskGlobal, headroomFrac, barFracFromTop, [0 0 1], pTxtGlobal, pTxtYOffsetFrac);
% Interaction bar (gray) shown BELOW the blue bar if significant
add_sig_bar_top_secondary(ax, timeVector, sigMaskInteract, barFracFromTop, 0.045, [0.50 0.50 0.50]);

% --- FIXED Y-SCALE FOR MANUSCRIPT FIGURES ---
ylim(ax, [-0.30 0.70]);

if zoomPlotsToStatsWin
    xlim(ax, statsWinMS);
end

legend(ax,[hGS hGD],'Location','southwest', FontSize=12);
grid(ax,'on');
box(ax,'off');
ax.TickDir = 'out';
ax.Layer   = 'top';
xline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');
%hold(ax,'off');
end

function add_sig_bar_top(ax, timeVector, sigMask, headroomFrac, barFracFromTop, rgb, pTxt, pTxtYOffsetFrac)
% add_sig_bar_top
% Draws a horizontal significance bar near the TOP of the axes (EEG style),
% and avoids covering the plotted curves by:
%   1) Adding headroom to ylim (upward expansion)
%   2) Placing the bar very near the top of the (expanded) axis range

if isempty(sigMask) || ~any(sigMask)
    return;
end

% ---------- Add headroom ----------
yl0 = ylim(ax);
yr0 = yl0(2) - yl0(1);
if yr0 <= 0
    yr0 = 1; % safety
end
yl = [yl0(1), yl0(2) + headroomFrac * yr0];
ylim(ax, yl);

yr = yl(2) - yl(1);

% ---------- Compute bar y-position (near top) ----------
yBar = yl(2) - barFracFromTop * yr;

% ---------- find contiguous TRUE segments ----------
sigMask = sigMask(:)';
d = diff([false sigMask false]);
starts = find(d==1);
ends   = find(d==-1)-1;

for i = 1:numel(starts)
    t1 = timeVector(starts(i));
    t2 = timeVector(ends(i));
    plot(ax, [t1 t2], [yBar yBar], '-', ...
        'Color', rgb, 'LineWidth', 4, 'HandleVisibility','off');
end

% ---------- p-text (optional) ----------
if nargin >= 7 && ~isempty(pTxt) && strlength(string(pTxt)) > 0
    xTxt = timeVector(max(1, starts(1))); % first cluster start
    yTxt = yBar + pTxtYOffsetFrac * yr;
    text(ax, xTxt, yTxt, char(pTxt), ...
        'Color', rgb, 'FontSize', 11, 'FontWeight','bold', ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
        'Interpreter','tex', 'Clipping','on', 'HandleVisibility','off');
end
end



%% ========================================================================
%  DUAL FIGURES (E1/E2 ONLY)  — keep these at the VERY END of the script
%  These functions ONLY affect the dual plots:
%    (1) plot_local_global_onefig_with_sig
%    (2) plot_type_E1_vs_E2_side_by_side
%
%  They place sig bars/labels relative to the SEM envelope (mu±sem),
%  and they add ONLY a small y-limit padding so bars/text never overwrite SEM.
% ========================================================================

function fig = plot_local_global_onefig_with_sig(timeVector, mu_local, sem_local, Ntr_local, ...
    mu_global, sem_global, Ntr_global, expLabel, Nsub, sigMaskLocal, sigMaskGlobal, sigMaskInteract, ...
    pTxtLocal, pTxtGlobal, statsWinMS, zoomPlotsToStatsWin, headroomFrac, barFracFromTop, pTxtYOffsetFrac, thisExp)
% 2-panel plot: Local (left) and Global (right), with TOP sig bars.
% Export size locked to match GRAND dual plot (3735 x 1247 @ 300 dpi).

% --- Fixed physical size so exportgraphics is consistent across monitors ---
FIG_W_IN = 3735/300;   % 12.45 in
FIG_H_IN = 1247/300;   % 4.1567 in

fig = figure('Color','w', 'Units','inches', 'Position',[0.5 0.5 FIG_W_IN FIG_H_IN], ...
    'InvertHardcopy','off');
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 FIG_W_IN FIG_H_IN]);
set(fig,'PaperSize',[FIG_W_IN FIG_H_IN]);

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
tl.Units = 'normalized';
tl.OuterPosition = [0 0 1 1];  % make tiles fill the figure canvas

% Colors (same as GRAND)
cStd = [0 0.45 0.85];
cDev = [0.90 0.45 0.10];

% ---------------- LEFT: Local ----------------
ax1 = nexttile; hold(ax1,'on');
set(ax1,'Color','w','XColor','k','YColor','k','FontSize',11);


        if strcmpi(thisExp,'E1')
            title(ax1, "Pupil size | Counting | Local effects", 'FontSize',12);
        else
            title(ax1, "Pupil size | Visual interference | Local effects", 'FontSize',12);
        end

            

xlabel(ax1,"Time (ms)", 'FontSize',12);
ylabel(ax1,"Pupil size (z-score)", 'FontSize',12);

upperLS = mu_local(1,:) + sem_local(1,:);
lowerLS = mu_local(1,:) - sem_local(1,:);
hf = fill(ax1,[timeVector fliplr(timeVector)],[upperLS fliplr(lowerLS)],cStd,'FaceAlpha',0.15,'EdgeColor','none');
set(hf,'HandleVisibility','off');
hLS = plot(ax1,timeVector,mu_local(1,:),'Color',cStd,'LineWidth',2,'DisplayName','Local Standard');

upperLD = mu_local(2,:) + sem_local(2,:);
lowerLD = mu_local(2,:) - sem_local(2,:);
hf = fill(ax1,[timeVector fliplr(timeVector)],[upperLD fliplr(lowerLD)],cDev,'FaceAlpha',0.15,'EdgeColor','none');
set(hf,'HandleVisibility','off');
hLD = plot(ax1,timeVector,mu_local(2,:),'Color',cDev,'LineWidth',2,'DisplayName','Local Deviant');

xline(ax1,0,'k-','LineWidth',1,'HandleVisibility','off');
grid(ax1,'on');
legend(ax1,[hLS hLD],'Location','southwest','FontSize',12);

if zoomPlotsToStatsWin
    xlim(ax1, statsWinMS);
end

% --- TOP sig bar placed above SEM (DATA COORDS) ---
lg_sigbar_top_from_sem(ax1, timeVector, sigMaskLocal, ...
    [upperLS; upperLD], [lowerLS; lowerLD], [0 0 1], pTxtLocal);
% Interaction bar (gray) shown BELOW the blue bar if significant
lg_sigbar_interaction_below_top_from_sem(ax1, timeVector, sigMaskInteract, ...
    [upperLS; upperLD], [lowerLS; lowerLD], [0.50 0.50 0.50]);

% --- FIXED Y-SCALE FOR MANUSCRIPT FIGURES ---
ylim(ax1, [-0.30 0.70]);

% ---------------- RIGHT: Global ----------------
ax2 = nexttile; hold(ax2,'on');
set(ax2,'Color','w','XColor','k','YColor','k','FontSize',11);


        if strcmpi(thisExp,'E1')
            title(ax2, "Pupil size | Counting | Global effects", 'FontSize',12);
        else
            title(ax2, "Pupil size | Visual interference | Global effects", 'FontSize',12);
        end



xlabel(ax2,"Time (ms)", 'FontSize',12);
ylabel(ax2,"Pupil size (z-score)", 'FontSize',12);

upperGS = mu_global(1,:) + sem_global(1,:);
lowerGS = mu_global(1,:) - sem_global(1,:);
hf = fill(ax2,[timeVector fliplr(timeVector)],[upperGS fliplr(lowerGS)],cStd,'FaceAlpha',0.15,'EdgeColor','none');
set(hf,'HandleVisibility','off');
hGS = plot(ax2,timeVector,mu_global(1,:),'Color',cStd,'LineWidth',2,'DisplayName','Global Standard');

upperGD = mu_global(2,:) + sem_global(2,:);
lowerGD = mu_global(2,:) - sem_global(2,:);
hf = fill(ax2,[timeVector fliplr(timeVector)],[upperGD fliplr(lowerGD)],cDev,'FaceAlpha',0.15,'EdgeColor','none');
set(hf,'HandleVisibility','off');
hGD = plot(ax2,timeVector,mu_global(2,:),'Color',cDev,'LineWidth',2,'DisplayName','Global Deviant');

xline(ax2,0,'k-','LineWidth',1,'HandleVisibility','off');
grid(ax2,'on');
legend(ax2,[hGS hGD],'Location','southwest','FontSize',12);

if zoomPlotsToStatsWin
    xlim(ax2, statsWinMS);
end

% --- TOP sig bar placed above SEM (DATA COORDS) ---
lg_sigbar_top_from_sem(ax2, timeVector, sigMaskGlobal, ...
    [upperGS; upperGD], [lowerGS; lowerGD], [0 0 1], pTxtGlobal);
% Interaction bar (gray) shown BELOW the blue bar if significant
lg_sigbar_interaction_below_top_from_sem(ax2, timeVector, sigMaskInteract, ...
    [upperGS; upperGD], [lowerGS; lowerGD], [0.50 0.50 0.50]);

% --- FIXED Y-SCALE FOR MANUSCRIPT FIGURES ---
ylim(ax2, [-0.30 0.70]);

end


function fig = plot_type_E1_vs_E2_side_by_side(COMPARE, timeVector, statsWinMS, zoomPlotsToStatsWin, ...
    headroomFrac, barFracFromTop, pTxtYOffsetFrac, pText) %#ok<INUSD>
% Two panels (E1 left, E2 right), with TOP(Local) and BOTTOM(Global) sig bars.
% Export size locked to match GRAND dual plot (3735 x 1247 @ 300 dpi).

required = {'E1','E2'};
for i = 1:numel(required)
    if ~isfield(COMPARE, required{i})
        error('COMPARE.%s is missing.', required{i});
    end
end

FIG_W_IN = 3735/300;
FIG_H_IN = 1247/300;

fig = figure('Color','w', 'Units','inches', 'Position',[0.5 0.5 FIG_W_IN FIG_H_IN], ...
    'InvertHardcopy','off');
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 FIG_W_IN FIG_H_IN]);
set(fig,'PaperSize',[FIG_W_IN FIG_H_IN]);

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
tl.Units = 'normalized';
tl.OuterPosition = [0 0 1 1];

colors = [ ...
    0 0.45 0.85;
    0 0.45 0.85;
    0.90 0.45 0.10;
    0.90 0.45 0.10];

lineStyles = {'--','-','--','-'};

    function [upperAll, lowerAll, lineH] = draw_panel(ax, expKey)
        hold(ax,'on');
        set(ax,'Color','w','XColor','k','YColor','k','FontSize',11);

        if strcmpi(expKey,'E1')
            title(ax, "Pupil size | Counting", 'FontSize',12);
        else
            title(ax, "Pupil size | Visual interference", 'FontSize',12);
        end

        xlabel(ax,"Time (ms)", 'FontSize',12);
        ylabel(ax,"Pupil size (z-score)", 'FontSize',12);

        n = numel(COMPARE.(expKey).typeOrder);
        upperAll = nan(n, numel(timeVector));
        lowerAll = nan(n, numel(timeVector));
        lineH    = gobjects(n,1);

        for k = 1:n
            mu  = COMPARE.(expKey).mu_type(k,:);
            sem = COMPARE.(expKey).sem_type(k,:);
            if all(isnan(mu)), continue; end

            upper = mu + sem;
            lower = mu - sem;

            upperAll(k,:) = upper;
            lowerAll(k,:) = lower;

            hFill = fill(ax, [timeVector fliplr(timeVector)], [upper fliplr(lower)], ...
                colors(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
            set(hFill,'HandleVisibility','off');

            lineH(k) = plot(ax, timeVector, mu, ...
                'Color', colors(k,:), 'LineStyle', lineStyles{k}, 'LineWidth', 2, ...
                'DisplayName', char(COMPARE.(expKey).typeOrder(k)));
        end

        xline(ax,0,'k-','LineWidth',1,'HandleVisibility','off');
        grid(ax,'on');

        if zoomPlotsToStatsWin
            xlim(ax, statsWinMS);
        end
    end

ax1 = nexttile; [up1, lo1, ~] = draw_panel(ax1,'E1');
ax2 = nexttile; [up2, lo2, lineH2] = draw_panel(ax2,'E2');

% Legend only on E2
legend(ax2, lineH2(isgraphics(lineH2)), "Location","northeast", "FontSize",10);


% =================== SIG BARS + LABELS (INSIDE FUNCTION) ===================

% ---- E1 labels ----
labLoc1 = "Local effect";
if pText && isfield(COMPARE.E1,'pTxtLocal') && strlength(string(COMPARE.E1.pTxtLocal))>0
    labLoc1 = "Local: " + string(COMPARE.E1.pTxtLocal);
end

labGlo1 = "Global effect";
if pText && isfield(COMPARE.E1,'pTxtGlobal') && strlength(string(COMPARE.E1.pTxtGlobal))>0
    labGlo1 = "Global: " + string(COMPARE.E1.pTxtGlobal);
end

% ---- E2 labels ----
labLoc2 = "Local effect";
if pText && isfield(COMPARE.E2,'pTxtLocal') && strlength(string(COMPARE.E2.pTxtLocal))>0
    labLoc2 = "Local: " + string(COMPARE.E2.pTxtLocal);
end

labGlo2 = "Global effect";
if pText && isfield(COMPARE.E2,'pTxtGlobal') && strlength(string(COMPARE.E2.pTxtGlobal))>0
    labGlo2 = "Global: " + string(COMPARE.E2.pTxtGlobal);
end

% ---- Draw bars/labels (DATA-COORDS, time-correct) ----
% TOP (Local): label placed BELOW the bar to avoid titles
lg_sigbar_top_from_sem(ax1, timeVector, COMPARE.E1.maskLocal,  up1, lo1, [0 0 1], labLoc1);
lg_sigbar_top_from_sem(ax2, timeVector, COMPARE.E2.maskLocal,  up2, lo2, [0 0 1], labLoc2);

% BOTTOM (Global): label placed ABOVE the bar to avoid x-label
lg_sigbar_bottom_from_sem(ax1, timeVector, COMPARE.E1.maskGlobal, up1, lo1, [0 0 1], labGlo1);
lg_sigbar_bottom_from_sem(ax2, timeVector, COMPARE.E2.maskGlobal, up2, lo2, [0 0 1], labGlo2);

% Interaction bar (gray) shown RIGHT ABOVE the blue global bar if significant
lg_sigbar_interaction_above_bottom_from_sem(ax1, timeVector, COMPARE.E1.maskInteract, up1, lo1, [0.50 0.50 0.50]);
lg_sigbar_interaction_above_bottom_from_sem(ax2, timeVector, COMPARE.E2.maskInteract, up2, lo2, [0.50 0.50 0.50]);

end


%% Helpers (dual plots only)

% function lg_dual_sigbar_top_from_sem(ax, timeVector, sigMask, upperCurves, lowerCurves, rgb, txt)
% % Draw TOP sig bar in DATA coords; y is placed just above SEM envelope.
% if isempty(sigMask) || ~any(sigMask), return; end
% 
% xl = xlim(ax);
% idxX = (timeVector >= xl(1)) & (timeVector <= xl(2));
% if ~any(idxX), idxX = true(size(timeVector)); end
% 
% u = upperCurves(:,idxX);
% l = lowerCurves(:,idxX);
% 
% yMax = max(u(:), [], 'omitnan');
% yMin = min(l(:), [], 'omitnan');
% if ~isfinite(yMax) || ~isfinite(yMin), return; end
% 
% yr = yMax - yMin; if yr <= 0, yr = 1; end
% barGap  = 0.04 * yr;
% textGap = 0.02 * yr;
% padTop  = 0.02 * yr;
% 
% yBar = yMax + barGap;
% yTxt = yBar + textGap;
% 
% yl = ylim(ax);
% if yTxt > yl(2)
%     ylim(ax, [yl(1), yTxt + padTop]);
% end
% 
% [starts, ends] = sig_segments(sigMask);
% for i = 1:numel(starts)
%     t1 = max(timeVector(starts(i)), xl(1));
%     t2 = min(timeVector(ends(i)),   xl(2));
%     if t2 <= t1, continue; end
%     plot(ax, [t1 t2], [yBar yBar], '-', 'Color', rgb, ...
%         'LineWidth', 4, 'HandleVisibility','off', 'Clipping','on');
% end
% 
% if nargin >= 7 && strlength(string(txt))>0
%     tTxt = min(max(timeVector(starts(1)), xl(1)), xl(2));
%     text(ax, tTxt, yTxt, char(txt), ...
%         'Color', rgb, 'FontSize', 11, 'FontWeight','bold', ...
%         'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
%         'Interpreter','tex', 'Clipping','on', 'HandleVisibility','off');
% end
% end


function lg_sigbar_top_from_sem(ax, timeVector, sigMask, upperCurves, lowerCurves, rgb, txt)
% Draw TOP sig bar in DATA coords, y-position set above SEM envelope.
% Text is placed BELOW the bar (prevents title overlap).

if isempty(sigMask) || ~any(sigMask), return; end

xl = xlim(ax);
idxX = (timeVector >= xl(1)) & (timeVector <= xl(2));
if ~any(idxX), idxX = true(size(timeVector)); end

u = upperCurves(:,idxX);
l = lowerCurves(:,idxX);

yMax = max(u(:), [], 'omitnan');
yMin = min(l(:), [], 'omitnan');
if ~isfinite(yMax) || ~isfinite(yMin), return; end

yr = yMax - yMin; 
if yr <= 0, yr = 1; end

barGap = 0.05 * yr;     % distance above SEM
txtGap = 0.03 * yr;     % distance below bar

yBar = yMax + barGap;
yTxt = yBar - txtGap;   % BELOW bar (key change)

% ensure axis includes the bar (not necessarily extra above it)
yl = ylim(ax);
if yBar > yl(2)
    ylim(ax, [yl(1), yBar + 0.03*yr]);
end

[starts, ends] = sig_segments(sigMask);
for i = 1:numel(starts)
    t1 = max(timeVector(starts(i)), xl(1));
    t2 = min(timeVector(ends(i)),   xl(2));
    if t2 <= t1, continue; end

    plot(ax, [t1 t2], [yBar yBar], '-', ...
        'Color', rgb, 'LineWidth', 4, ...
        'HandleVisibility','off', 'Clipping','on');
end

if nargin >= 7 && strlength(string(txt))>0
    tTxt = min(max(timeVector(starts(1)), xl(1)), xl(2));
    text(ax, tTxt, yTxt, char(txt), ...
        'Color', rgb, 'FontSize', 11, 'FontWeight','bold', ...
        'VerticalAlignment','top', 'HorizontalAlignment','left', ...
        'Interpreter','tex', 'Clipping','on', 'HandleVisibility','off');
end
end


function lg_sigbar_bottom_from_sem(ax, timeVector, sigMask, upperCurves, lowerCurves, rgb, txt)
% Draw BOTTOM sig bar in DATA coords, y-position set below SEM envelope.
% Text is placed ABOVE the bar (prevents x-label overlap).

if isempty(sigMask) || ~any(sigMask), return; end

xl = xlim(ax);
idxX = (timeVector >= xl(1)) & (timeVector <= xl(2));
if ~any(idxX), idxX = true(size(timeVector)); end

u = upperCurves(:,idxX);
l = lowerCurves(:,idxX);

yMax = max(u(:), [], 'omitnan');
yMin = min(l(:), [], 'omitnan');
if ~isfinite(yMax) || ~isfinite(yMin), return; end

yr = yMax - yMin; 
if yr <= 0, yr = 1; end

barGap = 0.06 * yr;     % distance below SEM
txtGap = 0.03 * yr;     % distance above bar

yBar = yMin - barGap;
yTxt = yBar + txtGap;   % ABOVE bar (key change)

yl = ylim(ax);
if yBar < yl(1)
    ylim(ax, [yBar - 0.03*yr, yl(2)]);
end

[starts, ends] = sig_segments(sigMask);
for i = 1:numel(starts)
    t1 = max(timeVector(starts(i)), xl(1));
    t2 = min(timeVector(ends(i)),   xl(2));
    if t2 <= t1, continue; end

    plot(ax, [t1 t2], [yBar yBar], '-', ...
        'Color', rgb, 'LineWidth', 4, ...
        'HandleVisibility','off', 'Clipping','on');
end

if nargin >= 7 && strlength(string(txt))>0
    tTxt = min(max(timeVector(starts(1)), xl(1)), xl(2));
    text(ax, tTxt, yTxt, char(txt), ...
        'Color', rgb, 'FontSize', 11, 'FontWeight','bold', ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
        'Interpreter','tex', 'Clipping','on', 'HandleVisibility','off');
end
end


function [starts, ends] = sig_segments(mask)
mask = mask(:)';                 
d = diff([false mask false]);
starts = find(d==1);
ends   = find(d==-1)-1;
end

function lg_print_png_fixed(fig, outPng, dpi)
% Always exports the FULL FIGURE at exact pixel size defined by PaperPosition/PaperSize.
% Use this ONLY for the dual figures to guarantee 3735x1247 at 300 dpi.

if nargin < 3 || isempty(dpi)
    dpi = 300;
end

% Ensure deterministic printing size
set(fig, 'InvertHardcopy','off');
set(fig, 'PaperPositionMode','manual');   % IMPORTANT (prevents MATLAB auto-scaling)

% Print respects PaperPosition/PaperSize (already set in your dual-plot functions)
print(fig, outPng, '-dpng', sprintf('-r%d', dpi));
end


%% ========================================================================
%  NEW 2-panel summary exports (PUPIL) in the SAME style as EEG exports
% ========================================================================
% Creates 5 fixed-size summary figures (3735x1247 px @300 dpi):
%   1) PUPIL_Bekinschtein_SSI_MATRIX_2panel_E1.png
%   2) PUPIL_Bekinschtein_SSI_MATRIX_2panel_E2.png
%   3) PUPIL_BinnedPvals_SubjectMatrix_2panel_E1.png
%   4) PUPIL_BinnedPvals_SubjectMatrix_2panel_E2.png
%   5) PUPIL_Bekinschtein_SSI_MEAN_INDEX_2panel_E1_E2.png
%
% Matches the EEG export style:
%   - same canvas
%   - same tiledlayout spacing
%   - same font sizing logic
%   - same matrix axis formatting
%   - same colorbar placement logic
%   - same "zoomed-out" look
% ========================================================================

function pupil_export_5_extra_zoomed_figures(outpath, grandOutFolder, SSI, cfgSSI, cfgBIN, excludeSubs, useField)

dpi = 300;
FIG_W_IN = 3735/dpi;
FIG_H_IN = 1247/dpi;

if ~exist(grandOutFolder,'dir'), error('grandOutFolder not found: %s', grandOutFolder); end
if ~exist(outpath,'dir'),        error('outpath not found: %s', outpath); end

%% ================================================================
% 1-2) SUBJECT-WISE SSI DARK MATRICES
% ================================================================
pb    = cfgSSI.pBins_matrix(:)';
nBins = numel(pb);
cmap  = pupil_pbin_cmap(cfgSSI, nBins+1);

% Post-0 selection
tms = SSI.E1.times(:)';
keepPost0 = true(size(tms));
if isfield(cfgSSI,'post0_only') && cfgSSI.post0_only
    keepPost0 = (tms >= 0);
end
tms2 = tms(keepPost0);

hideAbove = inf;
if isfield(cfgSSI,'displayHideBinsAboveP') && isfinite(cfgSSI.displayHideBinsAboveP)
    hideAbove = cfgSSI.displayHideBinsAboveP;
end

% ---------- E1 SSI ----------
if isfield(SSI,'E1') && ~isempty(SSI.E1.local)
    E1L = SSI.E1.local(:,keepPost0);
    E1G = SSI.E1.global(:,keepPost0);
    [E1L,E1G] = pupil_apply_ssi_hide(E1L,E1G,pb,hideAbove);

    fig = pupil_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    pupil_imagesc_matrix(ax1, tms2, E1L, ...
        'Pupil size | Counting', 'Local effects', cmap, nBins, 'Participants');

    ax2 = nexttile;
    pupil_imagesc_matrix(ax2, tms2, E1G, ...
        'Pupil size | Counting', 'Global effects', cmap, nBins, 'Participants');

    % match EEG SSI colorbar placement
    pupil_add_ssi_colorbar(ax2, pb, cmap, cfgSSI, [-0.5 nBins+0.5], [0.90 0.21 0.015 0.58]);

    outPng = fullfile(grandOutFolder,'PUPIL_Bekinschtein_SSI_MATRIX_2panel_E1.png');
    pupil_print_fixed(fig, outPng, dpi);
    close(fig);
end

% ---------- E2 SSI ----------
if isfield(SSI,'E2') && ~isempty(SSI.E2.local)
    E2L = SSI.E2.local(:,keepPost0);
    E2G = SSI.E2.global(:,keepPost0);
    [E2L,E2G] = pupil_apply_ssi_hide(E2L,E2G,pb,hideAbove);

    fig = pupil_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    pupil_imagesc_matrix(ax1, tms2, E2L, ...
        'Pupil size | Visual interference', 'Local effects', cmap, nBins, 'Participants');

    ax2 = nexttile;
    pupil_imagesc_matrix(ax2, tms2, E2G, ...
        'Pupil size | Visual interference', 'Global effects', cmap, nBins, 'Participants');

    % match EEG SSI colorbar placement
    pupil_add_ssi_colorbar(ax2, pb, cmap, cfgSSI, [-0.5 nBins+0.5], [0.90 0.21 0.015 0.58]);

    outPng = fullfile(grandOutFolder,'PUPIL_Bekinschtein_SSI_MATRIX_2panel_E2.png');
    pupil_print_fixed(fig, outPng, dpi);
    close(fig);
end

%% ================================================================
% 3-4) BINNED BLUE-WHITE MATRICES
% ================================================================
BIN    = pupil_build_binned_subject_pvals(outpath, cfgBIN, excludeSubs, useField);
cmapBW = pupil_cmap_bluewhite(256);

% ---------- E1 binned matrix ----------
if isfield(BIN,'E1') && ~isempty(BIN.E1.pLoc)
    fig = pupil_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    pupil_plot_binned_pvals(ax1, BIN.tCenters, BIN.E1.pLoc, BIN.E1.sigLoc, ...
        'Pupil size | Counting', 'Local effects', 'Participants', cmapBW, BIN.cMax);

    ax2 = nexttile;
    pupil_plot_binned_pvals(ax2, BIN.tCenters, BIN.E1.pGlo, BIN.E1.sigGlo, ...
        'Pupil size | Counting', 'Global effects', 'Participants', cmapBW, BIN.cMax);

    % match EEG binned-matrix colorbar placement
    cb = colorbar(ax2);
    cb.Position = [0.94 0.20 0.010 0.60];
    cb.Limits   = [0 BIN.cMax];
    cb.Label.String = 'p-values';
    cb.Label.Rotation = 90;
    cb.Label.Units = 'normalized';
    cb.Label.Position = [-1.35 0.50 0];
    cb.Label.HorizontalAlignment = 'center';
    cb.Label.VerticalAlignment   = 'middle';
    cb.FontSize = 12;

    outPng = fullfile(grandOutFolder,'PUPIL_BinnedPvals_SubjectMatrix_2panel_E1.png');
    pupil_print_fixed(fig, outPng, dpi);
    close(fig);
end

% ---------- E2 binned matrix ----------
if isfield(BIN,'E2') && ~isempty(BIN.E2.pLoc)
    fig = pupil_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    pupil_plot_binned_pvals(ax1, BIN.tCenters, BIN.E2.pLoc, BIN.E2.sigLoc, ...
        'Pupil size | Visual interference', 'Local effects', 'Participants', cmapBW, BIN.cMax);

    ax2 = nexttile;
    pupil_plot_binned_pvals(ax2, BIN.tCenters, BIN.E2.pGlo, BIN.E2.sigGlo, ...
        'Pupil size | Visual interference', 'Global effects', 'Participants', cmapBW, BIN.cMax);

    % match EEG binned-matrix colorbar placement
    cb = colorbar(ax2);
    cb.Position = [0.94 0.20 0.010 0.60];
    cb.Limits   = [0 BIN.cMax];
    cb.Label.String = 'p-values';
    cb.Label.Rotation = 90;
    cb.Label.Units = 'normalized';
    cb.Label.Position = [-1.35 0.50 0];
    cb.Label.HorizontalAlignment = 'center';
    cb.Label.VerticalAlignment   = 'middle';
    cb.FontSize = 12;

    outPng = fullfile(grandOutFolder,'PUPIL_BinnedPvals_SubjectMatrix_2panel_E2.png');
    pupil_print_fixed(fig, outPng, dpi);
    close(fig);
end

%% ================================================================
% 5) MEAN SSI INDEX LINE PLOT
% ================================================================
if isfield(SSI,'E1') && isfield(SSI,'E2') && ~isempty(SSI.E1.local) && ~isempty(SSI.E2.local)

    keep = keepPost0;
    tms2 = tms(keep);

    mE1L = mean(SSI.E1.local(:,keep),1,'omitnan');
    mE1G = mean(SSI.E1.global(:,keep),1,'omitnan');
    mE2L = mean(SSI.E2.local(:,keep),1,'omitnan');
    mE2G = mean(SSI.E2.global(:,keep),1,'omitnan');

    yMax = max([mE1L(:); mE1G(:); mE2L(:); mE2G(:)], [], 'omitnan');
    yMax = max(1, 1.05*yMax);

    fig = pupil_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.10 0.915 0.82];

    ax = nexttile;
    hold(ax,'on'); grid(ax,'on');
    set(ax,'FontSize',13,'LineWidth',0.8,'Box','off');
    plot(ax, tms2, mE1L, 'LineWidth',2);
    plot(ax, tms2, mE1G, 'LineWidth',2);
    title(ax,'Counting', 'FontSize',12, 'FontWeight','bold');
    xlabel(ax,'Time after 5th sound onset (ms)', 'FontSize',13);
    ylabel(ax,'Mean binned index', 'FontSize',13);
    % lgd = legend(ax,{'Local effect','Global effect'},'Location','northeast');
    % lgd.FontSize = 12;
    ylim(ax,[0 yMax]);

    ax = nexttile;
    hold(ax,'on'); grid(ax,'on');
    set(ax,'FontSize',13,'LineWidth',0.8,'Box','off');
    plot(ax, tms2, mE2L, 'LineWidth',2);
    plot(ax, tms2, mE2G, 'LineWidth',2);
    title(ax,'Visual interference', 'FontSize',12, 'FontWeight','bold');
    xlabel(ax,'Time after 5th sound onset (ms)', 'FontSize',13);
    ylabel(ax,'Mean binned index', 'FontSize',13);
    lgd = legend(ax,{'Local effect','Global effect'},'Location','northeast');
    lgd.FontSize = 12;
    ylim(ax,[0 yMax]);

    sgtitle(fig,'Mean individual statistical significance index | Pupil size','FontWeight','bold');

    outPng = fullfile(grandOutFolder,'PUPIL_Bekinschtein_SSI_MEAN_INDEX_2panel_E1_E2.png');
    pupil_print_fixed(fig, outPng, dpi);
    close(fig);
end

fprintf('[pupil_export_5_extra_zoomed_figures] DONE -> %s\n', grandOutFolder);
end

%% ======================================================================
% Helpers
% ======================================================================

function fig = pupil_fig_fixed(FIG_W_IN, FIG_H_IN)
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 FIG_W_IN FIG_H_IN], ...
    'InvertHardcopy','off');
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 FIG_W_IN FIG_H_IN]);
set(fig,'PaperSize',[FIG_W_IN FIG_H_IN]);
set(fig,'PaperPositionMode','manual');
end

function pupil_print_fixed(fig, outPng, dpi)
print(fig, outPng, '-dpng', sprintf('-r%d', dpi));
end

function cmap = pupil_pbin_cmap(cfgSSI, nColors)
if isfield(cfgSSI,'binColors') && ~isempty(cfgSSI.binColors) ...
        && size(cfgSSI.binColors,2)==3 && size(cfgSSI.binColors,1)==nColors
    cmap = cfgSSI.binColors;
else
    cmap = jet(max(nColors,2)); %#ok<JET>
end
end

function [A,B] = pupil_apply_ssi_hide(A,B,pBins,hideAboveP)
if ~isfinite(hideAboveP), return; end
hideLvl = find(pBins <= hideAboveP, 1, 'last');
if isempty(hideLvl), hideLvl = 0; end
A(A < hideLvl) = 0;
B(B < hideLvl) = 0;
end

function pupil_imagesc_matrix(ax, tms, M, titleTop, titleSmall, cmap, nBins, ylabStr)
imagesc(ax, tms, 1:size(M,1), M);

set(ax, ...
    'YDir','reverse', ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'Box','on', ...
    'Layer','top');

xlabel(ax,'Time after 5th sound onset (ms)', 'FontSize',13);
ylabel(ax,ylabStr, 'FontSize',13);

title(ax, sprintf('%s\n%s', titleTop, titleSmall), ...
    'FontSize',12, 'FontWeight','bold');

colormap(ax, cmap);
caxis(ax, [-0.5 nBins+0.5]);
grid(ax,'off');
end

function pupil_add_ssi_colorbar(ax, pBins, cmap, cfgSSI, clim, pos)
nBins = numel(pBins);
cb = colorbar(ax);
cb.Position = pos;
colormap(ax, cmap);
cb.Limits = clim;

ticks = 0:nBins;
cb.Ticks = ticks;

tickLabels = strings(size(ticks));
noneLabel = ">0.05";
if isfield(cfgSSI,'noneLabel') && strlength(string(cfgSSI.noneLabel))>0
    noneLabel = string(cfgSSI.noneLabel);
end
tickLabels(1) = noneLabel;
for k = 1:nBins
    tickLabels(k+1) = sprintf('<=%.4g', pBins(k));
end
cb.TickLabels = tickLabels;

lbl = 'p-values';
if isfield(cfgSSI,'binLabel') && strlength(string(cfgSSI.binLabel))>0
    lbl = char(cfgSSI.binLabel);
end
cb.Label.String = lbl;

% match EEG SSI label placement
cb.Label.Rotation = 90;
cb.Label.Units = 'normalized';
cb.Label.Position = [5.8 0.5 0];
cb.Label.HorizontalAlignment = 'center';
cb.Label.VerticalAlignment   = 'middle';
cb.FontSize = 12;
end

%% ======================================================================
% Binned p-values per subject
% ======================================================================

function BIN = pupil_build_binned_subject_pvals(outpath, cfgBIN, excludeSubs, useField)

BIN = struct();

binSize = cfgBIN.binSizeMS;
binStep = cfgBIN.binStepMS;
winMS   = cfgBIN.winMS;

tStarts  = winMS(1):binStep:(winMS(2)-binSize);
tCenters = tStarts + binSize/2;

BIN.tCenters = tCenters(:)';
BIN.cMax = 0.5;

for eidx = 1:2
    expKey = sprintf('E%d', eidx);

    subs = pupil_list_subjects(outpath, excludeSubs);
    pLoc = nan(numel(subs), numel(tCenters));
    pGlo = nan(numel(subs), numel(tCenters));

    keepRow = false(numel(subs),1);
    relTimes_ref = [];

    for i = 1:numel(subs)
        sID = subs(i);
        matFile = fullfile(outpath, char(sID), sprintf('%s_%s_pupil_preproc.mat', sID, expKey));
        if ~exist(matFile,'file'), continue; end

        S = load(matFile, 'dataTableClean','relTimes');
        if ~isfield(S,'dataTableClean') || isempty(S.dataTableClean) || height(S.dataTableClean)==0, continue; end

        relTimes = S.relTimes(:)';

        if isempty(relTimes_ref)
            relTimes_ref = relTimes;
        else
            if numel(relTimes)~=numel(relTimes_ref) || any(relTimes~=relTimes_ref)
                error('relTimes mismatch in %s (%s).', sID, expKey);
            end
        end

        T = S.dataTableClean;
        keepRow(i) = true;

        soundVec = string(T.sound);
        idxLS = (soundVec=="Local Standard");
        idxLD = (soundVec=="Local Deviant");

        typeVec = string(T.type);
        effectVec = strings(height(T),1);
        effectVec(startsWith(typeVec,"Frequent")) = "Global Standard";
        effectVec(startsWith(typeVec,"Rare"))     = "Global Deviant";
        idxGS = (effectVec=="Global Standard");
        idxGD = (effectVec=="Global Deviant");

        X = vertcat(T.(useField){:});

        XA = X(idxLS,:); XB = X(idxLD,:);
        XG = X(idxGS,:); XD = X(idxGD,:);

        for b = 1:numel(tCenters)
            wMask = (relTimes_ref >= tStarts(b)) & (relTimes_ref < (tStarts(b)+binSize));
            if ~any(wMask), continue; end

            a  = mean(XA(:,wMask), 2, 'omitnan');
            b2 = mean(XB(:,wMask), 2, 'omitnan');
            g  = mean(XG(:,wMask), 2, 'omitnan');
            d  = mean(XD(:,wMask), 2, 'omitnan');

            a  = a(~isnan(a)); 
            b2 = b2(~isnan(b2));
            g  = g(~isnan(g)); 
            d  = d(~isnan(d));

            if numel(a)>=2 && numel(b2)>=2
                [~, pLoc(i,b)] = ttest2(a,b2,'Vartype','unequal');
            end
            if numel(g)>=2 && numel(d)>=2
                [~, pGlo(i,b)] = ttest2(g,d,'Vartype','unequal');
            end
        end
    end

    subs = subs(keepRow);
    pLoc = pLoc(keepRow,:);
    pGlo = pGlo(keepRow,:);

    sigLoc = false(size(pLoc));
    sigGlo = false(size(pGlo));

    if isfield(cfgBIN,'doFWE') && cfgBIN.doFWE
        for r = 1:size(pLoc,1)
            sigLoc(r,:) = pupil_fwe_mask(pLoc(r,:), cfgBIN.alphaFWE, cfgBIN.fweMethod);
            sigGlo(r,:) = pupil_fwe_mask(pGlo(r,:), cfgBIN.alphaFWE, cfgBIN.fweMethod);
        end
    else
        sigLoc = (pLoc <= cfgBIN.alphaFWE);
        sigGlo = (pGlo <= cfgBIN.alphaFWE);
    end

    BIN.(expKey).subjects = subs;
    BIN.(expKey).pLoc     = pLoc;
    BIN.(expKey).pGlo     = pGlo;
    BIN.(expKey).sigLoc   = sigLoc;
    BIN.(expKey).sigGlo   = sigGlo;
end
end

function subs = pupil_list_subjects(outpath, excludeSubs)
d = dir(fullfile(outpath,'S*'));
d = d([d.isdir]);
names = string({d.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+', 'once'));
subs = names(keep);

if ~isempty(excludeSubs)
    subs = subs(~ismember(subs, excludeSubs));
end

nums = double(extractAfter(subs,"S"));
[~,ix] = sort(nums);
subs = subs(ix);
end

function sig = pupil_fwe_mask(pRow, alpha, method)
sig = false(size(pRow));
p = pRow(:);
ok = isfinite(p);
p2 = p(ok);

if isempty(p2), return; end

switch lower(string(method))
    case "bonferroni"
        sig(ok) = (p2 <= alpha/numel(p2));
    otherwise
        [ps,ix] = sort(p2,'ascend');
        m = numel(ps);
        thr = alpha ./ (m:-1:1)';
        k = find(ps > thr, 1, 'first');
        if isempty(k)
            keep = true(m,1);
        else
            keep = false(m,1);
            keep(1:k-1) = true;
        end
        sigOk = false(m,1);
        sigOk(ix) = keep;
        tmp = false(size(p));
        tmp(ok) = sigOk;
        sig = reshape(tmp, size(pRow));
end
end

function pupil_plot_binned_pvals(ax, tCenters, P, sigMask, titleTop, titleSmall, ylabStr, cmapBW, cMax)

Pdisp = P;
Pdisp(Pdisp > cMax) = cMax;
Pdisp(~isfinite(Pdisp)) = cMax;

imagesc(ax, tCenters, 1:size(Pdisp,1), Pdisp);

set(ax, ...
    'YDir','reverse', ...
    'FontSize',13, ...
    'LineWidth',0.8, ...
    'Box','on', ...
    'Layer','top');

xlabel(ax,'Time after 5th sound onset (ms)', 'FontSize',13);
ylabel(ax,ylabStr, 'FontSize',13);

title(ax, sprintf('%s\n%s', titleTop, titleSmall), ...
    'FontSize',12, 'FontWeight','bold');

colormap(ax, cmapBW);
caxis(ax, [0 cMax]);
grid(ax,'off');

[rr,cc] = find(sigMask);
hold(ax,'on');
plot(ax, tCenters(cc), rr, 'ro', 'MarkerSize',4, ...
    'MarkerFaceColor','r', 'LineWidth',0.5);
hold(ax,'off');
end

function cmap = pupil_cmap_bluewhite(n)
if nargin < 1 || isempty(n), n = 256; end

a = [247,251,255]/255; 
b = [222,235,247]/255; 
c = [198,219,239]/255;
d = [158,202,225]/255; 
e = [107,174,214]/255; 
f = [66,146,198]/255; 
g = [33,113,181]/255;

customColors = [g; f; e; d; c; b; a];
cmap = interp1(linspace(0,1,size(customColors,1)), ...
               customColors, linspace(0,1,n), 'pchip');
end

%% Interaction effects helpers

function add_sig_bar_top_secondary(ax, timeVector, sigMask, barFracFromTop, gapFrac, rgb)
% Draw a SECOND bar below the main top bar (used for interaction in 2-label plots)
if isempty(sigMask) || ~any(sigMask), return; end

yl = ylim(ax);
yr = yl(2) - yl(1);
if yr <= 0, yr = 1; end

% Main blue bar in add_sig_bar_top sits at:
%   yBlue = yl(2) - barFracFromTop * yr
yBlue = yl(2) - barFracFromTop * yr;

% Put gray bar a little below it
yBar = yBlue - gapFrac * yr;

[starts, ends] = sig_segments(sigMask);
for i = 1:numel(starts)
    t1 = timeVector(starts(i));
    t2 = timeVector(ends(i));
    plot(ax, [t1 t2], [yBar yBar], '-', ...
        'Color', rgb, 'LineWidth', 3.5, ...
        'HandleVisibility','off', 'Clipping','on');
end
end

function lg_sigbar_interaction_below_top_from_sem(ax, timeVector, sigMask, upperCurves, lowerCurves, rgb)
% Draw gray interaction bar BELOW the blue top bar in dual 2-label plots
if isempty(sigMask) || ~any(sigMask), return; end

xl = xlim(ax);
idxX = (timeVector >= xl(1)) & (timeVector <= xl(2));
if ~any(idxX), idxX = true(size(timeVector)); end

u = upperCurves(:,idxX);
l = lowerCurves(:,idxX);

yMax = max(u(:), [], 'omitnan');
yMin = min(l(:), [], 'omitnan');
if ~isfinite(yMax) || ~isfinite(yMin), return; end

yr = yMax - yMin;
if yr <= 0, yr = 1; end

% Blue top bar is at yMax + 0.05*yr in lg_sigbar_top_from_sem
% Put gray interaction bar below it, but still above SEM
yBar = yMax + 0.02 * yr;

yl = ylim(ax);
if yBar > yl(2)
    ylim(ax, [yl(1), yBar + 0.03*yr]);
end

[starts, ends] = sig_segments(sigMask);
for i = 1:numel(starts)
    t1 = max(timeVector(starts(i)), xl(1));
    t2 = min(timeVector(ends(i)),   xl(2));
    if t2 <= t1, continue; end

    plot(ax, [t1 t2], [yBar yBar], '-', ...
        'Color', rgb, 'LineWidth', 3.5, ...
        'HandleVisibility','off', 'Clipping','on');
end
end

function lg_sigbar_interaction_above_bottom_from_sem(ax, timeVector, sigMask, upperCurves, lowerCurves, rgb)
% Draw gray interaction bar ABOVE the blue bottom global bar in 4-label plots
if isempty(sigMask) || ~any(sigMask), return; end

xl = xlim(ax);
idxX = (timeVector >= xl(1)) & (timeVector <= xl(2));
if ~any(idxX), idxX = true(size(timeVector)); end

u = upperCurves(:,idxX);
l = lowerCurves(:,idxX);

yMax = max(u(:), [], 'omitnan');
yMin = min(l(:), [], 'omitnan');
if ~isfinite(yMax) || ~isfinite(yMin), return; end

yr = yMax - yMin;
if yr <= 0, yr = 1; end

% Blue bottom bar is at yMin - 0.06*yr in lg_sigbar_bottom_from_sem
% Put gray interaction bar just above it
yBar = yMin - 0.03 * yr;

yl = ylim(ax);
if yBar < yl(1)
    ylim(ax, [yBar - 0.03*yr, yl(2)]);
end

[starts, ends] = sig_segments(sigMask);
for i = 1:numel(starts)
    t1 = max(timeVector(starts(i)), xl(1));
    t2 = min(timeVector(ends(i)),   xl(2));
    if t2 <= t1, continue; end

    plot(ax, [t1 t2], [yBar yBar], '-', ...
        'Color', rgb, 'LineWidth', 3.5, ...
        'HandleVisibility','off', 'Clipping','on');
end
end

function S = summarize_interaction_for_paper(A_int, B_int, relTimes, statsWinMS, clustersInt, alphaFWE, expLabel)
% summarize_interaction_for_paper
% Creates a paper-style summary of the interaction:
%   A_int = FLD - FLS
%   B_int = RLD - RLS
%
% Returns:
%   S.table     -> one-row table
%   S.sentence  -> ready-to-copy sentence
%   S.t         -> timewise paired-t values
%   S.p         -> timewise paired p-values

relTimes = relTimes(:)';
idxWin   = (relTimes >= statsWinMS(1)) & (relTimes <= statsWinMS(2));
nTime    = numel(relTimes);

tVec = nan(1,nTime);
pVec = nan(1,nTime);

% ---------- Pointwise paired t-tests across subjects ----------
for tt = find(idxWin)
    a = A_int(:,tt);
    b = B_int(:,tt);

    keep = ~(isnan(a) | isnan(b));
    a = a(keep);
    b = b(keep);

    if numel(a) >= 2
        [~,p,~,st] = ttest(a,b);   % paired t-test
        tVec(tt) = st.tstat;
        pVec(tt) = p;
    end
end

% ---------- Pointwise summary ----------
ok = idxWin & isfinite(tVec) & isfinite(pVec);

if any(ok)
    idxOk = find(ok);
    [~,ix] = max(abs(tVec(ok)));   % use strongest |t|, report signed t there
    idxTmax   = idxOk(ix);
    tmax      = tVec(idxTmax);
    p_at_tmax = pVec(idxTmax);
    pmin      = min(pVec(ok));
    time_tmax = relTimes(idxTmax);
else
    tmax      = NaN;
    p_at_tmax = NaN;
    pmin      = NaN;
    time_tmax = NaN;
end

% ---------- Cluster summary ----------
nCandidateClusters = numel(clustersInt);

if nCandidateClusters > 0
    pFWE_all = [clustersInt.pFWE];
    mass_all = [clustersInt.clusterMass];

    [bestCluster_pFWE, iBest] = min(pFWE_all);
    bestCluster_mass = mass_all(iBest);

    startMS = relTimes(clustersInt(iBest).startIdx);
    endMS   = relTimes(clustersInt(iBest).endIdx);

    nSignificantClusters = sum(pFWE_all < alphaFWE);
    hasSignificantCluster = any(pFWE_all < alphaFWE);
else
    bestCluster_pFWE    = NaN;
    bestCluster_mass    = NaN;
    startMS             = NaN;
    endMS               = NaN;
    nSignificantClusters = 0;
    hasSignificantCluster = false;
end

% ---------- Ready-to-write sentence ----------
if hasSignificantCluster
    sentence = sprintf(['Interaction effect: significant cluster ' ...
        'pFWE = %.3f, tmax = %.2f, p at tmax = %.3f, startMS = %.0f, endMS = %.0f'], ...
        bestCluster_pFWE, tmax, p_at_tmax, startMS, endMS);

elseif nCandidateClusters > 0
    sentence = sprintf(['No interaction effect: no significant cluster ' ...
        'pFWE = %.3f, tmax = %.2f, p at tmax = %.3f, startMS = %.0f, endMS = %.0f'], ...
        bestCluster_pFWE, tmax, p_at_tmax, startMS, endMS);

else
    sentence = sprintf(['No interaction effect: no candidate cluster, ' ...
        'tmax = %.2f, p at tmax = %.3f'], ...
        tmax, p_at_tmax);
end

% ---------- One-row table ----------
S.table = table( ...
    string(expLabel), ...
    hasSignificantCluster, ...
    nCandidateClusters, ...
    nSignificantClusters, ...
    bestCluster_pFWE, ...
    bestCluster_mass, ...
    startMS, ...
    endMS, ...
    time_tmax, ...
    tmax, ...
    p_at_tmax, ...
    pmin, ...
    string(sentence), ...
    'VariableNames', { ...
        'Experiment', ...
        'HasSignificantCluster', ...
        'NumCandidateClusters', ...
        'NumSignificantClusters', ...
        'BestCluster_pFWE', ...
        'BestCluster_Mass', ...
        'StartMS', ...
        'EndMS', ...
        'Time_tmax_MS', ...
        'tmax', ...
        'p_at_tmax', ...
        'pmin_uncorrected', ...
        'ReportSentence'});

S.sentence = sentence;
S.t        = tVec;
S.p        = pVec;
end

function Tinteraction = export_interaction_summary_table(STATS, grandOutFolder)
% export_interaction_summary_table
% Collects E1/E2 interaction summary rows into one single table,
% saves them in GRAND_AVG, and also saves a TXT file with ready-made sentences.

rows = {};

expList = fieldnames(STATS);
for i = 1:numel(expList)
    expKey = expList{i};
    if isfield(STATS.(expKey), 'interactionSummary') && isfield(STATS.(expKey).interactionSummary, 'table')
        rows{end+1,1} = STATS.(expKey).interactionSummary.table; %#ok<AGROW>
    end
end

if isempty(rows)
    warning('No interactionSummary found in STATS. No table was exported.');
    Tinteraction = table();
    return;
end

Tinteraction = vertcat(rows{:});

% Save CSV
outCsv = fullfile(grandOutFolder, 'STATS_PUPIL_interaction_summary.csv');
writetable(Tinteraction, outCsv);

% Save XLSX (if available)
outXlsx = fullfile(grandOutFolder, 'STATS_PUPIL_interaction_summary.xlsx');
try
    writetable(Tinteraction, outXlsx);
catch
    warning('Could not save XLSX file. CSV was saved successfully.');
end

% Save TXT with copy-paste sentences
outTxt = fullfile(grandOutFolder, 'STATS_PUPIL_interaction_summary_sentences.txt');
fid = fopen(outTxt, 'w');
if fid ~= -1
    for r = 1:height(Tinteraction)
        fprintf(fid, '%s: %s\n', Tinteraction.Experiment(r), Tinteraction.ReportSentence(r));
    end
    fclose(fid);
else
    warning('Could not create TXT summary file.');
end

fprintf('[INTERACTION] Saved CSV:  %s\n', outCsv);
fprintf('[INTERACTION] Saved TXT:  %s\n', outTxt);
if exist(outXlsx, 'file')
    fprintf('[INTERACTION] Saved XLSX: %s\n', outXlsx);
end
end
