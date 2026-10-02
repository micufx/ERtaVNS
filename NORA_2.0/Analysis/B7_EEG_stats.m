clc; clear; close all;

% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% EEG_bekinschtein_stats_only_PAPER.m  (+ optional Monte-Carlo permutation validation)
% Stats-only pipeline (does NOT touch your grand-average building code).
% Loads single-subject ALL sets:
%   out\Sxx\auditory_all_LG_Sxx_E1.set  and  ..._E2.set
%
% Reconstructs the 4 LG bins and derives (FACTORIAL, equal-weight):
%   Local:  Standard = mean(FLS, RLS), Deviant = mean(FLD, RLD)
%   Global: Standard = mean(FLS, FLD), Deviant = mean(RLS, RLD)
%
% Bekinschtein 2009-style stats logic:
%   GROUP: paired t-tests across subjects (Dev vs Std) at each chan/time
%          + extent (triple) criterion:
%             p<=0.01 for >=minConsec consecutive samples on >=minElec electrodes
%   INDIVIDUAL: unpaired t-tests across trials within subject (Dev vs Std)
%          + extent criterion:
%             p<=0.01 for >=minConsec consecutive samples on >=minElec electrodes
%
% Optional Monte-Carlo (paired swap within-subject):
%   Validates whether "the extent criterion happens anywhere" is rare under H0.
%   If pMC >= cfg.mc.alpha: you may choose to reject the GROUP mask (default here: YES)
%
% IMPORTANT: This script uses the PAPER "count electrodes" rule.
% NO neighbor/cluster contiguity is used anywhere.
%
% Miguel Contreras-Altamirano, 2026

%% Paths
mainpath  = 'L:\Downloads\eeglab2023.0\';
outroot   = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';
grandFold = fullfile(outroot, 'GRAND_AVG');
outpath = 'L:\Downloads\scripts_Vagus\matlab\testing\out\';
rasFile = 'L:\Downloads\scripts_Vagus\matlab\testing\out\RAS_Local_Global.csv';

if ~exist(grandFold,'dir'), mkdir(grandFold); end

if exist('eeglab','file') ~= 2
    addpath(mainpath);
end

%% EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; %#ok<ASGLU>

%% Experiments
expList = {'E1','E2'};

%% Your 4 LG labels (as in EEG_grand_avg)
targets = ["Frequent Local Standard","Rare Local Standard", ...
    "Frequent Local Deviant","Rare Local Deviant"]; % [FLS RLS FLD RLD]

%% OPTIONAL SKIP LISTS
skipSubjects_all = {"S21","S12","S23","S02","S29"};
skipSubjects_byExp = struct();
skipSubjects_byExp.E1 = {};
skipSubjects_byExp.E2 = {};

%% Config (paper-consistent)
cfg = struct();

cfg.tlim_plot  = [-2500 700];
cfg.tlim_stats = [0 700];

cfg.showFromMs = 0; % display bars only after 0 ms

cfg.erpChans = {'Fz','FC1','Cz','Pz'};
cfg.ylim = [-8 8];

% Snap times from paper Fig.2
cfg.snapTimes_ms = [100 132 164 196 228 260 292 324 356 388 420 452 484];

% --- PAPER THRESHOLDS ---

% ERP group average
cfg.alpha_grp     = 0.05; % Bekinschtein uses p<=0.01 for the sample-wise criterion
cfg.minConsec_grp = 20;   % (paper uses 10 consecutive samples = 40ms)   at 500Hz == 20 samples
cfg.minElec_grp   = 5;    % adapted to 32ch (paper used 10 electrodes) --> 15.6% of total electrodes

% ERP individual level
cfg.alpha_ind     = 0.05;  % Bekinschtein uses p<=0.01 for the sample-wise criterion/individual extent-forming threshold
cfg.minConsec_ind = 20;    % (paper uses 5 consecutive = 20ms)  at 500Hz == 10 samples
cfg.minElec_ind   = 5;     % adapted to 32ch (paper used 10 electrodes) --> 15.6% of total electrodes

% ===================== DISPLAY BINS (match PUPIL) =====================
cfg.pBins_matrix = [0.05 0.01 0.005 0.001 0.0005 0.0001];  % EDIT
cfg.noneLabel    = ">0.05";
cfg.binLabel     = "p-values";

% ---- COLORS (match PUPIL: level 0 = black, then dark-blue -> red) ----
a        = [0 0 0];               % level 0 (>0.05)
blueDark = [33,113,181]/255;
redDark  = [178,24,43]/255;

nBins  = numel(cfg.pBins_matrix); % 6
gradBR = interp1([0 1], [blueDark; redDark], linspace(0,1,nBins), 'pchip'); % 6 colors

% binColors must be (nBins+1) x 3 : level0 + levels1..nBins
cfg.binColors = [a; gradBR];

cfg.displayHideBinsAboveP = inf;

cfg.caxis_erpTopo = [-2 2];

cfg.topoGrid     = 120;
cfg.topoContours = 6;
cfg.topoShading  = 'interp';
cfg.topoShading_stats = 'interp';
cfg.statsTopoMapLimits = [0 1];

cfg.exportDPI = 300;

cfg.doStratifiedBalancing = true;
cfg.balanceSeed = 7;

cfg.indTtest = "student"; % "student" (paper-safe) or "welch"

% Bar display (display only)
cfg.barMinConsec   = 0;
cfg.barYoffset_uV  = 1.2;
cfg.barMaxGap      = 0;

% Topos per condition
cfg.caxis_erpTopo_cond = [-2 2];  % example


%% ===================== SIGNIFICANCE BAR TEXT (NEW) =====================
% If false: blue bars still drawn, but NO p-text is printed.
% (Same behavior as your PUPIL script)
cfg.pText = false;   % true / false

%% ===================== MULTIPLE COMPARISON CORRECTION (MCC) =====================
% Choose ONE:
%   'none'            : no correction (paper-like)
%   'fdr_bh'          : Benjamini-Hochberg FDR control
%   'holm'            : Holm-Bonferroni (FWER)
%   'bonferroni'      : Bonferroni (FWER)
%   'fwe_permutation' : permutation FWE using max-|t| (GROUP paired maps only)
%
% Notes:
% - MCC is applied to p-maps BEFORE the triple-criterion.
% - Monte-Carlo EXTENT validation (cfg.mc.*) remains separate and runs AFTER masks.
% - If you select 'fwe_permutation', it will only be applied to GROUP paired maps.
cfg.mcc = struct();
cfg.mcc.method = 'fwe_permutation';   % <-- CHANGE HERE
cfg.mcc.q      = 0.05;     % for 'fdr_bh'
cfg.mcc.alpha  = 0.05;     % for 'holm'/'bonferroni' (not strictly needed; p-map stays p-adj)
cfg.mcc.fwe = struct();    % for 'fwe_permutation'
cfg.mcc.fwe.nPerm   = 10000;
cfg.mcc.fwe.seed    = 777;
cfg.mcc.fwe.verbose = true;
cfg.mcc.fwe.progressStep = 200;

% Where to apply MCC:
%   'none'       : never apply MCC (even if method ~= 'none')
%   'individual' : only individual (trial-level) p-maps
%   'group'      : only group (subject-level) p-maps
%   'both'       : apply to both individual + group
cfg.mcc.applyTo = 'group';   % <-- CHANGE THIS (recommended default: 'group')


%% ===================== BINNED EEG STATS (NEW) =====================
% ===== EEG BINNED MATRICES (PUPIL-STYLE) =====
% ===================== EEG BINNED MATRIX EXPORT (PUPIL-STYLE) =====================
doEEGBinned = true;

cfgBIN_EEG = struct();
cfgBIN_EEG.winMS     = cfg.tlim_stats;     % e.g. [0 700]
cfgBIN_EEG.binSizeMS = 50;                 % 50 or 100
cfgBIN_EEG.binStepMS = 25;                 % 50 (overlap if <binSize)

cfgBIN_EEG.alpha_ind   = cfg.alpha_ind;    % e.g. 0.01
cfgBIN_EEG.minElec_ind = cfg.minElec_ind;  % e.g. 5

cfgBIN_EEG.doFWE     = true;
cfgBIN_EEG.fweMethod = "holm";             % "holm" or "bonferroni"
cfgBIN_EEG.alphaFWE  = 0.05;

cfgBIN_EEG.climVals  = [0 0.05];           % good for EEG; change if you want
cfgBIN_EEG.indTtest  = cfg.indTtest;       % "student" or "welch"
cfgBIN_EEG.doBalance = true;
cfgBIN_EEG.balanceSeed = cfg.balanceSeed;

% optional:
cfgBIN_EEG.maskHeatmap = false;            % false = pupil-like full heatmap
cfgBIN_EEG.minConsecBins_ind = 20;          % 1 = no extent; >1 adds extent in bins

if doEEGBinned
    [outSubPng, outElecPng] = eeg_export_binned_trial_matrices(outroot, grandFold, cfgBIN_EEG, skipSubjects_all, targets);
    fprintf('Saved EEG binned matrices:\n  %s\n  %s\n', outSubPng, outElecPng);
end


%% ===================== MONTE-CARLO PERMUTATION VALIDATION =====================
cfg.mc = struct();
cfg.mc.doMC      = false;   % set false to see raw paper-like extent mask only
cfg.mc.nPerm     = 5000;
cfg.mc.alpha     = 0.05;
cfg.mc.seed      = 123;
cfg.mc.verbose   = true;
cfg.mc.progressStep = 100;

% If TRUE: reject (zero) group masks when MC does not validate them
cfg.mc.zeroMaskIfNotValidated = true;
% ================================================================================

%% Find participants
dSubs = dir(outroot);
isSub = [dSubs.isdir] & ~ismember({dSubs.name},{'.','..','GRAND_AVG'});
subNames = {dSubs(isSub).name};
subNames = subNames(~cellfun('isempty', regexp(subNames, '^S\d+', 'once')));

participant_list = {};
for si = 1:numel(subNames)
    p = subNames{si};
    fAllE1 = fullfile(outroot, p, sprintf('auditory_all_LG_%s_E1.set', p));
    if exist(fAllE1,'file')
        participant_list{end+1} = p; %#ok<SAGROW>
    end
end
participant_list = sort(participant_list);
fprintf('Participants found (>=E1 ALL): %d\n', numel(participant_list));
if isempty(participant_list), error('No participants found with expected ALL sets.'); end

%% Normalize skip lists
skip_all_norm = unique(lower(strtrim(string(skipSubjects_all))));

%% Hold SSI + GROUP p-maps
SSI = struct();
GROUPPMAP = struct();

for eidx = 1:numel(expList)
    thisExp = expList{eidx};
    fprintf('\n=== STATS for %s ===\n', thisExp);

    skipNow = skip_all_norm;
    if isfield(skipSubjects_byExp, thisExp)
        skipNow = [skipNow, lower(strtrim(string(skipSubjects_byExp.(thisExp))))];
    end
    skipNow = unique(skipNow);

    if ~isempty(skipNow)
        fprintf('Skip list for %s: %s\n', thisExp, strjoin(cellstr(skipNow), ', '));
    end

    %% Load template
    EEGtmpl = [];
    for si = 1:numel(participant_list)
        p = participant_list{si};
        fAll = fullfile(outroot, p, sprintf('auditory_all_LG_%s_%s.set', p, thisExp));
        if exist(fAll,'file')
            EEGtmpl = pop_loadset(fAll);
            EEGtmpl = eeg_checkset(EEGtmpl);
            break;
        end
    end
    if isempty(EEGtmpl)
        warning('No data found for %s. Skipping.', thisExp);
        continue;
    end

    times = double(EEGtmpl.times(:))';
    tidx_plot  = find(times >= cfg.tlim_plot(1)  & times <= cfg.tlim_plot(2));
    tidx_stats = find(times >= cfg.tlim_stats(1) & times <= cfg.tlim_stats(2));

    tms_plot  = times(tidx_plot);
    tms_stats = times(tidx_stats);

    %% Preallocate subject-level ERPs: [sub x chan x time]
    subIDs = strings(0,1);
    Lstd = []; Ldev = [];
    Gstd = []; Gdev = [];
    SSI_local  = [];
    SSI_global = [];

    nTrialLocalStd  = 0; nTrialLocalDev  = 0;
    nTrialGlobalStd = 0; nTrialGlobalDev = 0;

    % --- Mazancieux 4-condition storage (per-subject ERPs) ---
    FLSsub = [];  % Frequent Local Standard
    RLSsub = [];  % Rare Local Standard
    FLDsub = [];  % Frequent Local Deviant
    RLDsub = [];  % Rare Local Deviant


    %% Subject loop
    for si = 1:numel(participant_list)

        participant = participant_list{si};
        rng(cfg.balanceSeed + si);

        if ismember(lower(string(participant)), skipNow)
            fprintf('   [%s %s] In skip list -> skipping subject.\n', participant, thisExp);
            continue;
        end

        fAll = fullfile(outroot, participant, sprintf('auditory_all_LG_%s_%s.set', participant, thisExp));
        if ~exist(fAll,'file')
            fprintf('   [%s %s] missing ALL set. Skipping.\n', participant, thisExp);
            continue;
        end

        EEGs = [];
        while isempty(EEGs)
            try
                EEGs = pop_loadset(fAll);
                EEGs = eeg_checkset(EEGs);
            catch
                warning('pop_loadset failed for %s. Retrying once...', fAll);
                pause(2);
                EEGs = [];
            end
        end

        if ~(isfield(EEGs,'epoch') && isfield(EEGs.epoch,'eventlg_type'))
            fprintf('   [%s %s] missing EEG.epoch.eventlg_type. Skipping.\n', participant, thisExp);
            continue;
        end

        EEGs = reorder_to_template(EEGs, {EEGtmpl.chanlocs.labels});

        if numel(EEGs.times) ~= numel(EEGtmpl.times) || any(abs(double(EEGs.times(:)) - double(EEGtmpl.times(:))) > 1e-6)
            fprintf('   [%s %s] times mismatch. Skipping.\n', participant, thisExp);
            continue;
        end

        % Build 4 condition masks per trial
        M = false(4, EEGs.trials);
        for ep = 1:EEGs.trials
            lablist = to_strlist(EEGs.epoch(ep).eventlg_type);
            for k = 1:4
                M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
            end
        end

        idxFLS = find(M(1,:));
        idxRLS = find(M(2,:));
        idxFLD = find(M(3,:));
        idxRLD = find(M(4,:));

        if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
            fprintf('   [%s %s] missing one LG bin (FLS/RLS/FLD/RLD). Skipping.\n', participant, thisExp);
            continue;
        end

        % SUBJECT ERPs (FACTORIAL, equal weighting)
        ERP_FLS = mean(double(EEGs.data(:,:,idxFLS)), 3);
        ERP_RLS = mean(double(EEGs.data(:,:,idxRLS)), 3);
        ERP_FLD = mean(double(EEGs.data(:,:,idxFLD)), 3);
        ERP_RLD = mean(double(EEGs.data(:,:,idxRLD)), 3);

        lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3);
        ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3);
        gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3);
        gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3);

        subIDs(end+1,1) = string(participant); %#ok<SAGROW>
        Lstd(end+1,:,:) = lstd; %#ok<SAGROW>
        Ldev(end+1,:,:) = ldev; %#ok<SAGROW>
        Gstd(end+1,:,:) = gstd; %#ok<SAGROW>
        Gdev(end+1,:,:) = gdev; %#ok<SAGROW>

        nTrialLocalStd  = nTrialLocalStd  + numel(idxFLS) + numel(idxRLS);
        nTrialLocalDev  = nTrialLocalDev  + numel(idxFLD) + numel(idxRLD);
        nTrialGlobalStd = nTrialGlobalStd + numel(idxFLS) + numel(idxFLD);
        nTrialGlobalDev = nTrialGlobalDev + numel(idxRLS) + numel(idxRLD);


        % --- keep elementary condition ERPs for Mazancieux-style plotting ---
        FLSsub(end+1,:,:) = ERP_FLS; %#ok<SAGROW>
        RLSsub(end+1,:,:) = ERP_RLS; %#ok<SAGROW>
        FLDsub(end+1,:,:) = ERP_FLD; %#ok<SAGROW>
        RLDsub(end+1,:,:) = ERP_RLD; %#ok<SAGROW>


        % INDIVIDUAL trial-level stats
        if cfg.doStratifiedBalancing
            rng(cfg.balanceSeed);
            [devLoc, stdLoc] = build_balanced_trials(EEGs.data, idxFLD, idxRLD, idxFLS, idxRLS);
            [devGlo, stdGlo] = build_balanced_trials(EEGs.data, idxRLS, idxRLD, idxFLS, idxFLD);
        else
            devLoc = double(EEGs.data(:,:, [idxFLD(:); idxRLD(:)] ));
            stdLoc = double(EEGs.data(:,:, [idxFLS(:); idxRLS(:)] ));
            devGlo = double(EEGs.data(:,:, [idxRLS(:); idxRLD(:)] ));
            stdGlo = double(EEGs.data(:,:, [idxFLS(:); idxFLD(:)] ));
        end

        fprintf('   [%s %s] computing pLoc...\n', participant, thisExp);
        pLoc = unpaired_pmap(devLoc(:,tidx_stats,:), stdLoc(:,tidx_stats,:), cfg.indTtest);

        fprintf('   [%s %s] computing pGlo...\n', participant, thisExp);
        pGlo = unpaired_pmap(devGlo(:,tidx_stats,:), stdGlo(:,tidx_stats,:), cfg.indTtest);


        % --- MCC on individual p-maps (optional; controlled by cfg.mcc.applyTo) ---
        mccInfoIndLoc = struct('method','skipped');
        mccInfoIndGlo = struct('method','skipped');

        if should_apply_mcc(cfg, 'individual')
            [pLoc, mccInfoIndLoc] = apply_mcc_individual(devLoc(:,tidx_stats,:), stdLoc(:,tidx_stats,:), pLoc, cfg, 'IND_LOCAL', cfg.indTtest);
            [pGlo, mccInfoIndGlo] = apply_mcc_individual(devGlo(:,tidx_stats,:), stdGlo(:,tidx_stats,:), pGlo, cfg, 'IND_GLOBAL', cfg.indTtest);
        end


        maskLoc_ind = triple_criterion_time_mask(pLoc, cfg.alpha_ind, cfg.minElec_ind, cfg.minConsec_ind);
        maskGlo_ind = triple_criterion_time_mask(pGlo, cfg.alpha_ind, cfg.minElec_ind, cfg.minConsec_ind);

        % SSI_local(end+1,:)  = p_to_binned_index_timecourse(pLoc, maskLoc_ind, cfg.pBins_matrix);
        % SSI_global(end+1,:) = p_to_binned_index_timecourse(pGlo, maskGlo_ind, cfg.pBins_matrix);

        SSI_local(end+1,:) = p_to_binned_index_timecourse(pLoc, maskLoc_ind, cfg.pBins_matrix, cfg.minElec_ind);
        SSI_global(end+1,:) = p_to_binned_index_timecourse(pGlo, maskGlo_ind, cfg.pBins_matrix, cfg.minElec_ind);


        fprintf('   [%s %s] done p-maps.\n', participant, thisExp);
    end

    nSub = size(Lstd,1);
    if nSub < 2
        warning('Not enough subjects for %s (n=%d). Skipping.', thisExp, nSub);
        continue;
    end

    fprintf('-> %s: Nsub=%d | Local Std=%d Dev=%d | Global Std=%d Dev=%d\n', ...
        thisExp, nSub, nTrialLocalStd, nTrialLocalDev, nTrialGlobalStd, nTrialGlobalDev);

    %% GROUP stats (paired across subjects)
    pLoc_grp = paired_pmap(Ldev(:,:,tidx_stats), Lstd(:,:,tidx_stats));
    pGlo_grp = paired_pmap(Gdev(:,:,tidx_stats), Gstd(:,:,tidx_stats));

    % --- MCC on GROUP p-maps (paired across subjects) ---
    % This does NOT touch your Monte-Carlo EXTENT validation (cfg.mc.*). Separate steps.
    % If method == 'fwe_permutation', we recompute an FWE-corrected p-map using max-|t| sign-flips.

    mccInfoLoc = struct('method','skipped');
    mccInfoGlo = struct('method','skipped');

    if should_apply_mcc(cfg, 'group')
        [pLoc_grp, mccInfoLoc] = apply_mcc_group(Ldev(:,:,tidx_stats), Lstd(:,:,tidx_stats), pLoc_grp, cfg, 'GROUP_LOCAL');
        [pGlo_grp, mccInfoGlo] = apply_mcc_group(Gdev(:,:,tidx_stats), Gstd(:,:,tidx_stats), pGlo_grp, cfg, 'GROUP_GLOBAL');
    end


    % Optional: print what was applied
    if isfield(mccInfoLoc,'method')
        fprintf('[%s] MCC GROUP Local: %s\n', thisExp, mccInfoLoc.method);
    end
    if isfield(mccInfoGlo,'method')
        fprintf('[%s] MCC GROUP Global: %s\n', thisExp, mccInfoGlo.method);
    end


    % Raw (paper-like) extent masks BEFORE MC
    maskLoc_raw = triple_criterion_time_mask(pLoc_grp, cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp);
    maskGlo_raw = triple_criterion_time_mask(pGlo_grp, cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp);

    maskLoc_grp = maskLoc_raw;
    maskGlo_grp = maskGlo_raw;

    % ===================== MONTE-CARLO VALIDATION (GROUP) =====================
    pMC_Loc = NaN; nHit_Loc = NaN;
    pMC_Glo = NaN; nHit_Glo = NaN;

    if cfg.mc.doMC
        fprintf('\n[%s] Monte-Carlo validation (nPerm=%d)...\n', thisExp, cfg.mc.nPerm);

        % Only validate if there is something to validate
        obsHasLoc = any(maskLoc_raw);
        obsHasGlo = any(maskGlo_raw);

        if obsHasLoc
            [pMC_Loc, nHit_Loc] = mc_validate_triplecriterion_group( ...
                Ldev(:,:,tidx_stats), Lstd(:,:,tidx_stats), ...
                cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp, ...
                cfg.mc.nPerm, cfg.mc.seed + 10*eidx + 1, cfg.mc.verbose, cfg.mc.progressStep);
        else
            pMC_Loc = 1; nHit_Loc = cfg.mc.nPerm;
        end

        if obsHasGlo
            [pMC_Glo, nHit_Glo] = mc_validate_triplecriterion_group( ...
                Gdev(:,:,tidx_stats), Gstd(:,:,tidx_stats), ...
                cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp, ...
                cfg.mc.nPerm, cfg.mc.seed + 10*eidx + 2, cfg.mc.verbose, cfg.mc.progressStep);
        else
            pMC_Glo = 1; nHit_Glo = cfg.mc.nPerm;
        end

        fprintf('[%s] MC: Local pMC=%.4f (hits=%d/%d) | Global pMC=%.4f (hits=%d/%d)\n', ...
            thisExp, pMC_Loc, nHit_Loc, cfg.mc.nPerm, pMC_Glo, nHit_Glo, cfg.mc.nPerm);

        if cfg.mc.zeroMaskIfNotValidated
            if pMC_Loc >= cfg.mc.alpha
                fprintf('[%s] MC: Local NOT validated (pMC>=%.2f) -> zeroing GROUP Local mask.\n', thisExp, cfg.mc.alpha);
                maskLoc_grp(:) = false;
            end
            if pMC_Glo >= cfg.mc.alpha
                fprintf('[%s] MC: Global NOT validated (pMC>=%.2f) -> zeroing GROUP Global mask.\n', thisExp, cfg.mc.alpha);
                maskGlo_grp(:) = false;
            end
        end
    end


    % ---- store MCC metadata ----
    GROUPPMAP.(thisExp).mcc_method = char(lower(string(cfg.mcc.method)));
    if isfield(cfg.mcc,'q'), GROUPPMAP.(thisExp).mcc_q = cfg.mcc.q; else, GROUPPMAP.(thisExp).mcc_q = NaN; end
    if isfield(cfg.mcc,'fwe') && isfield(cfg.mcc.fwe,'nPerm'), GROUPPMAP.(thisExp).mcc_fwe_nPerm = cfg.mcc.fwe.nPerm; else, GROUPPMAP.(thisExp).mcc_fwe_nPerm = NaN; end
    if isfield(cfg.mcc,'fwe') && isfield(cfg.mcc.fwe,'seed'),  GROUPPMAP.(thisExp).mcc_fwe_seed  = cfg.mcc.fwe.seed;  else, GROUPPMAP.(thisExp).mcc_fwe_seed  = NaN; end

    if isfield(cfg.mcc,'applyTo')
        GROUPPMAP.(thisExp).mcc_applyTo = char(lower(string(cfg.mcc.applyTo)));
    else
        GROUPPMAP.(thisExp).mcc_applyTo = 'both'; % default assumption
    end

    % ---- quick flags: any significant timepoint? ----
    GROUPPMAP.(thisExp).obsHasLocal_raw  = any(maskLoc_raw);
    GROUPPMAP.(thisExp).obsHasGlobal_raw = any(maskGlo_raw);
    GROUPPMAP.(thisExp).obsHasLocal_postMC  = any(maskLoc_grp);
    GROUPPMAP.(thisExp).obsHasGlobal_postMC = any(maskGlo_grp);


    % Store group p-maps + both raw and post-MC masks
    GROUPPMAP.(thisExp).times = tms_stats;
    GROUPPMAP.(thisExp).pLoc  = pLoc_grp;
    GROUPPMAP.(thisExp).pGlo  = pGlo_grp;
    GROUPPMAP.(thisExp).mLoc_raw = maskLoc_raw;
    GROUPPMAP.(thisExp).mGlo_raw = maskGlo_raw;
    GROUPPMAP.(thisExp).mLoc  = maskLoc_grp;
    GROUPPMAP.(thisExp).mGlo  = maskGlo_grp;
    GROUPPMAP.(thisExp).pMC_local  = pMC_Loc;
    GROUPPMAP.(thisExp).pMC_global = pMC_Glo;
    GROUPPMAP.(thisExp).nHit_local  = nHit_Loc;
    GROUPPMAP.(thisExp).nHit_global = nHit_Glo;
    % ============================================================================


    %% Interaction effects
    INT.(thisExp) = eeg_export_interaction_effects( ...
        EEGtmpl, FLSsub, RLSsub, FLDsub, RLDsub, cfg, thisExp, grandFold);


    %% ERP figures
    outLocal = fullfile(grandFold, sprintf('Bekinschtein_STATS_%s_ERP_LOCAL.png', thisExp));
    outGlobal= fullfile(grandFold, sprintf('Bekinschtein_STATS_%s_ERP_GLOBAL.png', thisExp));

    [tf, loc] = ismember(tidx_plot, tidx_stats);


    % --- interaction masks mapped to plot window (for grey bars) ---
    % --- TRUE OMNIBUS interaction mask mapped to plot window (single grey bar source) ---
    pIntOmni_plot = ones(size(INT.(thisExp).grey_omnibus.p,1), numel(tidx_plot));
    mIntOmni_plot = false(1, numel(tidx_plot));

    pIntOmni_plot(:, tf) = INT.(thisExp).grey_omnibus.p(:, loc(tf));
    mIntOmni_plot(1, tf) = INT.(thisExp).grey_omnibus.mask(loc(tf));



    pLoc_plot = ones(size(pLoc_grp,1), numel(tidx_plot));
    mLoc_plot = false(1, numel(tidx_plot));
    pLoc_plot(:, tf) = pLoc_grp(:, loc(tf));
    mLoc_plot(1, tf) = maskLoc_grp(loc(tf));

    plot_erps_2x2(EEGtmpl, tms_plot, ...
        Lstd(:,:,tidx_plot), Ldev(:,:,tidx_plot), ...
        pLoc_plot, mLoc_plot, cfg, thisExp, 'LOCAL', nSub, nTrialLocalStd, nTrialLocalDev, outLocal, ...
        pIntOmni_plot, mIntOmni_plot);


    pGlo_plot = ones(size(pGlo_grp,1), numel(tidx_plot));
    mGlo_plot = false(1, numel(tidx_plot));
    pGlo_plot(:, tf) = pGlo_grp(:, loc(tf));
    mGlo_plot(1, tf) = maskGlo_grp(loc(tf));

    plot_erps_2x2(EEGtmpl, tms_plot, ...
        Gstd(:,:,tidx_plot), Gdev(:,:,tidx_plot), ...
        pGlo_plot, mGlo_plot, cfg, thisExp, 'GLOBAL', nSub, nTrialGlobalStd, nTrialGlobalDev, outGlobal, ...
        pIntOmni_plot, mIntOmni_plot);


    % --- Mazancieux-style 4-condition figure (Local top bar + Global bottom bar) ---
    outMzCombined = fullfile(grandFold, sprintf('Mazancieux_STATS_%s_ERP4_LOCAL_GLOBAL.png', thisExp));

    plot_mz4_2x2_with_dual_stats(EEGtmpl, tms_plot, ...
        FLSsub(:,:,tidx_plot), RLSsub(:,:,tidx_plot), FLDsub(:,:,tidx_plot), RLDsub(:,:,tidx_plot), ...
        pLoc_plot, mLoc_plot, pGlo_plot, mGlo_plot, ...
        cfg, thisExp, nSub, outMzCombined, ...
        pIntOmni_plot, mIntOmni_plot);



    %% EXTRA pupil-style EEG exports (2 channels per figure; fixed 3735x1247) =====

    eeg_export_pupilstyle_2panel_extras(EEGtmpl, tms_plot, ...
        Lstd(:,:,tidx_plot), Ldev(:,:,tidx_plot), ...
        Gstd(:,:,tidx_plot), Gdev(:,:,tidx_plot), ...
        pLoc_plot, mLoc_plot, pGlo_plot, mGlo_plot, ...
        FLSsub(:,:,tidx_plot), RLSsub(:,:,tidx_plot), FLDsub(:,:,tidx_plot), RLDsub(:,:,tidx_plot), ...
        cfg, thisExp, nSub, ...
        nTrialLocalStd, nTrialLocalDev, nTrialGlobalStd, nTrialGlobalDev, ...
        grandFold, ...
        pIntOmni_plot, mIntOmni_plot);


    %% ERP topo snaps

    % Topos effects (difference between conditions)
    diffLoc = squeeze(mean(Ldev(:,:,tidx_stats)-Lstd(:,:,tidx_stats),1));
    diffGlo = squeeze(mean(Gdev(:,:,tidx_stats)-Gstd(:,:,tidx_stats),1));

    outTopoLocERP = fullfile(grandFold, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_LOCAL.png', thisExp));
    outTopoGloERP = fullfile(grandFold, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_GLOBAL.png', thisExp));
    plot_topo_snaps_erp(EEGtmpl, tms_stats, diffLoc, cfg, sprintf('%s | ERP scalp maps (Dev-Std) | LOCAL', thisExp), outTopoLocERP);
    plot_topo_snaps_erp(EEGtmpl, tms_stats, diffGlo, cfg, sprintf('%s | ERP scalp maps (Dev-Std) | GLOBAL', thisExp), outTopoGloERP);

    % Topos per condition
    eeg_export_topo_snaps_conditions(EEGtmpl, tms_stats, ...
        Lstd(:,:,tidx_stats), Ldev(:,:,tidx_stats), ...
        Gstd(:,:,tidx_stats), Gdev(:,:,tidx_stats), ...
        cfg, thisExp, grandFold);

    %% STATS topo snaps
    outTopoLocST = fullfile(grandFold, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_LOCAL.png', thisExp));
    outTopoGloST = fullfile(grandFold, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_GLOBAL.png', thisExp));
    plot_topo_snaps_stats(EEGtmpl, tms_stats, pLoc_grp, maskLoc_grp, cfg, sprintf('%s | STATS scalp maps (p-values) | LOCAL', thisExp), outTopoLocST);
    plot_topo_snaps_stats(EEGtmpl, tms_stats, pGlo_grp, maskGlo_grp, cfg, sprintf('%s | STATS scalp maps (p-values) | GLOBAL', thisExp), outTopoGloST);


    %% EXTRA FULL-SCREEN STATS TOPO SNAPS | transparent background
    % This only adds extra figures. It does NOT change the statistics.

    outTopoLocST_full = fullfile(grandFold, ...
        sprintf('Bekinschtein_STATS_%s_TOPO_STATS_LOCAL_FULLSCREEN_NOBG.png', thisExp));

    outTopoGloST_full = fullfile(grandFold, ...
        sprintf('Bekinschtein_STATS_%s_TOPO_STATS_GLOBAL_FULLSCREEN_NOBG.png', thisExp));

    plot_topo_snaps_stats_fullscreen_nobg( ...
        EEGtmpl, tms_stats, pLoc_grp, maskLoc_grp, cfg, ...
        sprintf('%s | STATS scalp maps (p-values) | LOCAL | red dots: p<=0.05 (only if triple-criterion time)', thisExp), ...
        outTopoLocST_full);

    plot_topo_snaps_stats_fullscreen_nobg( ...
        EEGtmpl, tms_stats, pGlo_grp, maskGlo_grp, cfg, ...
        sprintf('%s | STATS scalp maps (p-values) | GLOBAL | red dots: p<=0.05 (only if triple-criterion time)', thisExp), ...
        outTopoGloST_full);



    %% Store SSI
    SSI.(thisExp).subIDs = subIDs;
    SSI.(thisExp).times  = tms_stats;
    SSI.(thisExp).local  = SSI_local;
    SSI.(thisExp).global = SSI_global;

end


%% Regularity Awareness Score
OUT = eeg_export_total_scatter_4ch(outpath, grandFold, GROUPPMAP, rasFile, skipSubjects_all);


%% SUBJECT-WISE MATRIX + mean
if isfield(SSI,'E1') && isfield(SSI,'E2')
    outMat = fullfile(grandFold, 'Bekinschtein_STATS_SSI_MATRIX.png');
    plot_subjectwise_matrix_post0(SSI, cfg, outMat);

    outMean = fullfile(grandFold, 'Bekinschtein_STATS_mean_individual_SSI_E1_vs_E2.png');
    plot_mean_subjectwise_index_post0(SSI, cfg, outMean);
end

%% ELECTRODE-WISE GROUP matrix
if isfield(GROUPPMAP,'E1') && isfield(GROUPPMAP,'E2')
    outGrp = fullfile(grandFold, 'Bekinschtein_STATS_GROUP_PVAL_MATRIX_2x2.png');
    plot_group_pval_matrix_2x2(GROUPPMAP, cfg, outGrp);
end


%% ===== NEW: export ALL 2-panel summary figures (9 total) =====
%eeg2_export_all_2panel_summary_figures(grandFold, SSI, GROUPPMAP, cfg, cfgBIN_EEG);

eeg_export_9_extra_zoomed_figures(outroot, grandFold, SSI, GROUPPMAP, cfg, cfgBIN_EEG, skipSubjects_all, targets);

%% MONTE-CARLO SUMMARY TABLE (MATLAB table + Excel copy)
if isfield(cfg,'mc')
    outXLS = fullfile(grandFold, 'Bekinschtein_STATS_MC_summary.xlsx');
    outCSV = fullfile(grandFold, 'Bekinschtein_STATS_MC_summary.csv');
    outMAT = fullfile(grandFold, 'Bekinschtein_STATS_MC_summary.mat');

    exps = intersect(fieldnames(GROUPPMAP), {'E1','E2'}, 'stable');

    Exp = strings(0,1);
    nSub = NaN(0,1);

    alpha_grp = NaN(0,1);
    minElec_grp = NaN(0,1);
    minConsec_grp = NaN(0,1);

    doMC = false(0,1);
    nPerm = NaN(0,1);
    alphaMC = NaN(0,1);

    obsHasLocal  = false(0,1);
    obsHasGlobal = false(0,1);

    pMC_local  = NaN(0,1);
    pMC_global = NaN(0,1);

    keepLocal  = false(0,1);
    keepGlobal = false(0,1);

    for i = 1:numel(exps)
        E = exps{i};
        Exp(end+1,1) = string(E);

        if isfield(SSI, E) && isfield(SSI.(E),'subIDs')
            nSub(end+1,1) = numel(SSI.(E).subIDs);
        else
            nSub(end+1,1) = NaN;
        end

        alpha_grp(end+1,1)     = cfg.alpha_grp;
        minElec_grp(end+1,1)   = cfg.minElec_grp;
        minConsec_grp(end+1,1) = cfg.minConsec_grp;

        doMC(end+1,1)     = cfg.mc.doMC;
        nPerm(end+1,1)    = cfg.mc.nPerm;
        alphaMC(end+1,1)  = cfg.mc.alpha;

        obsHasLocal(end+1,1)  = isfield(GROUPPMAP.(E),'mLoc_raw') && any(GROUPPMAP.(E).mLoc_raw);
        obsHasGlobal(end+1,1) = isfield(GROUPPMAP.(E),'mGlo_raw') && any(GROUPPMAP.(E).mGlo_raw);

        pMC_local(end+1,1)  = getfield_default(GROUPPMAP.(E),'pMC_local', NaN);
        pMC_global(end+1,1) = getfield_default(GROUPPMAP.(E),'pMC_global', NaN);

        keepLocal(end+1,1)  = isfield(GROUPPMAP.(E),'mLoc') && any(GROUPPMAP.(E).mLoc);
        keepGlobal(end+1,1) = isfield(GROUPPMAP.(E),'mGlo') && any(GROUPPMAP.(E).mGlo);
    end

    Tmc = table(Exp, nSub, alpha_grp, minElec_grp, minConsec_grp, ...
        doMC, nPerm, alphaMC, obsHasLocal, obsHasGlobal, pMC_local, pMC_global, keepLocal, keepGlobal);

    save(outMAT, 'Tmc');
    writetable(Tmc, outXLS, 'FileType','spreadsheet');
    writetable(Tmc, outCSV);

    disp('Monte-Carlo summary table (Tmc):');
    disp(Tmc);
end


%% ===================== RUN SUMMARY (MCC + EXTENT + MC) =====================
% Prints a compact report:
% - MCC method used
% - extent rule parameters
% - whether anything was significant before MC (raw extent)
% - MC p-values + whether masks survived after MC
% - sample counts if available

try
    print_run_summary(GROUPPMAP, SSI, cfg);
catch
    warning('Could not print run summary.');
end
% ============================================================================


%% Slope Pz ERP
OUT_EEG_RAINBOX_PZ = eeg_export_rainbox_compare_E1E2_Pz( ...
    outroot, grandFold, GROUPPMAP, skipSubjects_all, targets, 'Pz');

SUPP_EEG = eeg_export_supplementary_tables_all( ...
    outroot, grandFold, GROUPPMAP, SSI, cfg, ...
    skipSubjects_all, targets, ["Pz","Fz"], OUT_EEG_RAINBOX_PZ, INT);


disp('Stats pipeline done.');

%% FINAL RESULTS

% Window times summary effects

fprintf('\nN subjects:\n');
fprintf('E1: %d\n', numel(SSI.E1.subIDs));
fprintf('E2: %d\n', numel(SSI.E2.subIDs));

for E = ["E1","E2"]
    fprintf('\n=== %s ===\n', E);

    t = GROUPPMAP.(E).times(:)';

    % Local
    m = GROUPPMAP.(E).mLoc(:)';
    d = diff([false m false]);
    rs = find(d==1);
    re = find(d==-1)-1;

    if isempty(rs)
        fprintf('Local: no significant window\n');
    else
        for k = 1:numel(rs)
            fprintf('Local window %d: %.0f to %.0f ms\n', ...
                k, t(rs(k)), t(re(k)));
        end
    end

    % Global
    m = GROUPPMAP.(E).mGlo(:)';
    d = diff([false m false]);
    rs = find(d==1);
    re = find(d==-1)-1;

    if isempty(rs)
        fprintf('Global: no significant window\n');
    else
        for k = 1:numel(rs)
            fprintf('Global window %d: %.0f to %.0f ms\n', ...
                k, t(rs(k)), t(re(k)));
        end
    end

    if isfield(GROUPPMAP.(E),'pMC_local')
        fprintf('pMC local: %.4f\n', GROUPPMAP.(E).pMC_local);
    end
    if isfield(GROUPPMAP.(E),'pMC_global')
        fprintf('pMC global: %.4f\n', GROUPPMAP.(E).pMC_global);
    end
end


print_run_summary(GROUPPMAP, SSI, cfg)


%%
% Channel specific effects summary

roiLabels = {'Pz','Fz'};
expKeys   = {'E1','E2'};
effects   = {'Loc','Glo'};

for e = 1:numel(expKeys)
    expKey = expKeys{e};

    for ef = 1:numel(effects)
        thisEff = effects{ef};

        if strcmpi(thisEff,'Loc')
            thisMask = GROUPPMAP.(expKey).mLoc;
            thisP    = GROUPPMAP.(expKey).pLoc;
            effName  = 'Local effect';
        else
            thisMask = GROUPPMAP.(expKey).mGlo;
            thisP    = GROUPPMAP.(expKey).pGlo;
            effName  = 'Global effect';
        end

        fprintf('\n=== %s | %s ===\n', expKey, effName);

        for i = 1:numel(roiLabels)
            ch = find(strcmpi({EEGtmpl.chanlocs.labels}, roiLabels{i}), 1);

            if isempty(ch)
                fprintf('%s not found.\n', roiLabels{i});
                continue;
            end

            sigIdx = thisMask & (thisP(ch,:) <= cfg.alpha_grp);

            if any(sigIdx)
                pvals = thisP(ch,sigIdx);
                tvals = GROUPPMAP.(expKey).times(sigIdx);

                [pmin, ix] = min(pvals);
                tmin = tvals(ix);

                fprintf('%s | %s | min pFWE = %.8f at %g ms\n', ...
                    roiLabels{i}, effName, pmin, tmin);
            else
                fprintf('%s | %s | no significant samples in ROI window\n', ...
                    roiLabels{i}, effName);
            end
        end
    end
end

%%

% Interaction effect results summary
for expName = ["E1","E2"]
    fprintf('\n=== %s INTERACTION ===\n', expName);

    G = INT.(expName).grey_omnibus;

    if isfield(G,'times')
        t = G.times(:)';
    elseif isfield(INT.(expName),'times')
        t = INT.(expName).times(:)';
    elseif isfield(GROUPPMAP.(expName),'times')
        t = GROUPPMAP.(expName).times(:)';
    else
        error('%s interaction: no time vector found.', expName);
    end

    p = G.p;
    m = logical(G.mask(:)');

    if size(p,2) ~= numel(t)
        error('%s interaction: size mismatch between p-map columns and time vector.', expName);
    end

    if ~any(m)
        fprintf('No significant interaction window\n');
        continue;
    end

    d  = diff([false m false]);
    rs = find(d == 1);
    re = find(d == -1) - 1;

    for k = 1:numel(rs)
        winIdx = rs(k):re(k);

        fprintf('Interaction window %d: %.0f to %.0f ms\n', ...
            k, t(rs(k)), t(re(k)));

        pWin = p(:, winIdx);
        [pmin, linIdx] = min(pWin(:));
        [chIdx, tIdxLocal] = ind2sub(size(pWin), linIdx);
        tIdx = winIdx(tIdxLocal);

        fprintf('  min pFWE = %.8f at %s, %.0f ms\n', ...
            pmin, EEGtmpl.chanlocs(chIdx).labels, t(tIdx));
    end
end

%% Individual stats summary
fprintf('\n=== INDIVIDUAL-LEVEL EEG SUMMARY ===\n');

for expName = ["E1","E2"]
    t = SSI.(expName).times(:)';
    keep = t >= 0 & t <= 700;
    t2 = t(keep);
    dt = median(diff(t2));

    fprintf('\n%s\n', expName);

    for effName = ["local","global"]
        X = SSI.(expName).(effName)(:, keep);   % subjects x time, binned index 0..6

        nSub = size(X,1);
        hasAny = any(X > 0, 2);
        nAny = sum(hasAny);
        pctAny = 100 * nAny / nSub;

        dur_ms = sum(X > 0, 2) * dt;

        mSSI = mean(X, 1, 'omitnan');
        [peakSSI, ixPeak] = max(mSSI);
        peakTime = t2(ixPeak);

        fprintf('%s | subjects with any effect: %d/%d (%.1f%%)\n', ...
            char(effName), nAny, nSub, pctAny);
        fprintf('%s | peak mean SSI: %.3f at %.0f ms\n', ...
            char(effName), peakSSI, peakTime);
        fprintf('%s | median duration among all subjects: %.0f ms\n', ...
            char(effName), median(dur_ms, 'omitnan'));
        fprintf('%s | median duration among subjects with effect: %.0f ms\n', ...
            char(effName), median(dur_ms(hasAny), 'omitnan'));
    end
end


%% Interaction effects summary 
%% ============================================================
% CHANNEL-SPECIFIC INTERACTION BARS: Pz and Fz
% Prints exactly the samples represented by the grey bars
% =============================================================

roiLabels = ["Pz","Fz"];
chLabels  = string({EEGtmpl.chanlocs.labels});

fprintf('\n============================================================\n');
fprintf('CHANNEL-SPECIFIC LOCAL-BY-GLOBAL INTERACTION BARS\n');
fprintf('============================================================\n');

for E = ["E1","E2"]

    tms  = double(GROUPPMAP.(E).times(:)');
    pInt = double(INT.(E).grey_omnibus.p);
    mInt = logical(INT.(E).grey_omnibus.mask(:)');

    fprintf('\n%s\n', E);

    for r = 1:numel(roiLabels)

        roi = roiLabels(r);
        ch  = find(strcmpi(chLabels, roi), 1);

        if isempty(ch)
            fprintf('%s: channel not found\n', roi);
            continue;
        end

        % This reproduces the grey-bar logic used in the figure
        showMask = mInt & ...
                   (tms >= cfg.showFromMs) & ...
                   (pInt(ch,:) <= cfg.alpha_grp);

        showMask = eeg_bridge_small_gaps(showMask, cfg.barMaxGap);
        showMask = eeg_keep_min_run_length(showMask, cfg.barMinConsec);

        if ~any(showMask)
            fprintf('%s: no corrected interaction samples\n', roi);
            continue;
        end

        d  = diff([false showMask false]);
        rs = find(d == 1);
        re = find(d == -1) - 1;

        for k = 1:numel(rs)

            idx = rs(k):re(k);

            [pMin, localIdx] = min(pInt(ch,idx));
            minIdx = idx(localIdx);

            fprintf(['%s | interval %d: %.0f to %.0f ms | ' ...
                     'min pFWE = %.8f at %.0f ms\n'], ...
                roi, k, tms(rs(k)), tms(re(k)), ...
                pMin, tms(minIdx));
        end
    end
end


%% ================================================================
% QUANTITATIVE ROI SUMMARY FOR MAIN EEG RESULTS TEXT
% Adds window-based and peak-sample descriptive statistics for Pz/Fz
% using the SAME subject files and SAME local/global reconstruction.
%
% IMPORTANT:
% - pFWE is still your PRIMARY corrected result.
% - The t-tests below are FOLLOW-UP ROI summaries on window-averaged
%   amplitudes (or on the peak sample), useful for reporting magnitudes.
% ================================================================

fprintf('\n============================================================\n');
fprintf('QUANTITATIVE ROI SUMMARY FOR MAIN EEG RESULTS\n');
fprintf('============================================================\n');

roiLabels = {'Pz','Fz'};
expKeys   = {'E1','E2'};

for e = 1:numel(expKeys)
    expKey = expKeys{e};

    if ~isfield(SSI, expKey) || ~isfield(SSI.(expKey), 'subIDs') || isempty(SSI.(expKey).subIDs)
        fprintf('\n%s: no subjects found in SSI.\n', expKey);
        continue;
    end

    % --- load one template file for channel lookup ---
    subjList = cellstr(SSI.(expKey).subIDs);
    EEGtmp = [];
    for s = 1:numel(subjList)
        fTmp = fullfile(outroot, subjList{s}, sprintf('auditory_all_LG_%s_%s.set', subjList{s}, expKey));
        if exist(fTmp,'file')
            EEGtmp = pop_loadset(fTmp);
            EEGtmp = eeg_checkset(EEGtmp);
            break;
        end
    end
    if isempty(EEGtmp)
        fprintf('\n%s: no template EEG file found.\n', expKey);
        continue;
    end

    tStats = GROUPPMAP.(expKey).times(:)';

    for eff = ["Loc","Glo"]

        if eff == "Loc"
            effLong = 'Local';
            thisMask = GROUPPMAP.(expKey).mLoc(:)';
            thisP    = GROUPPMAP.(expKey).pLoc;
        else
            effLong = 'Global';
            thisMask = GROUPPMAP.(expKey).mGlo(:)';
            thisP    = GROUPPMAP.(expKey).pGlo;
        end

        fprintf('\n------------------------------------------------------------\n');
        fprintf('%s | %s effect\n', expKey, effLong);
        fprintf('------------------------------------------------------------\n');

        d = diff([false thisMask false]);
        rs = find(d == 1);
        re = find(d == -1) - 1;

        if isempty(rs)
            fprintf('No significant scalp-level window for %s.\n', effLong);
            continue;
        end

        for w = 1:numel(rs)

            tStart = tStats(rs(w));
            tEnd   = tStats(re(w));

            fprintf('\nWindow %d: %.0f to %.0f ms\n', w, tStart, tEnd);

            for r = 1:numel(roiLabels)

                roi = roiLabels{r};
                chTmp = find(strcmpi({EEGtmp.chanlocs.labels}, roi), 1);

                if isempty(chTmp)
                    fprintf('\n%s | channel not found in template.\n', roi);
                    continue;
                end

                % ROI-specific significant samples INSIDE the scalp-level window
                roiSigIdx = false(size(thisMask));
                roiSigIdx(rs(w):re(w)) = thisMask(rs(w):re(w)) & (thisP(chTmp, rs(w):re(w)) <= cfg.alpha_grp);

                hasROIsig = any(roiSigIdx);

                if hasROIsig
                    pvalsROI = thisP(chTmp, roiSigIdx);
                    tvalsROI = tStats(roiSigIdx);
                    [pminROI, ixMin] = min(pvalsROI);
                    tPeakROI = tvalsROI(ixMin);
                else
                    pminROI = NaN;
                    tPeakROI = NaN;
                end

                % subject-level storage
                yStdWin  = nan(numel(subjList),1);
                yDevWin  = nan(numel(subjList),1);
                yStdPeak = nan(numel(subjList),1);
                yDevPeak = nan(numel(subjList),1);

                for s = 1:numel(subjList)

                    subj = subjList{s};
                    fAll = fullfile(outroot, subj, sprintf('auditory_all_LG_%s_%s.set', subj, expKey));
                    if ~exist(fAll,'file')
                        continue;
                    end

                    EEGs = pop_loadset(fAll);
                    EEGs = eeg_checkset(EEGs);

                    if ~(isfield(EEGs,'epoch') && isfield(EEGs.epoch,'eventlg_type'))
                        continue;
                    end

                    ch = find(strcmpi({EEGs.chanlocs.labels}, roi), 1);
                    if isempty(ch)
                        continue;
                    end

                    M = false(4, EEGs.trials);
                    for ep = 1:EEGs.trials
                        lablist = to_strlist(EEGs.epoch(ep).eventlg_type);
                        for k = 1:4
                            M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
                        end
                    end

                    idxFLS = find(M(1,:));
                    idxRLS = find(M(2,:));
                    idxFLD = find(M(3,:));
                    idxRLD = find(M(4,:));

                    if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
                        continue;
                    end

                    FLS = squeeze(mean(double(EEGs.data(ch,:,idxFLS)),3)); FLS = FLS(:)';
                    RLS = squeeze(mean(double(EEGs.data(ch,:,idxRLS)),3)); RLS = RLS(:)';
                    FLD = squeeze(mean(double(EEGs.data(ch,:,idxFLD)),3)); FLD = FLD(:)';
                    RLD = squeeze(mean(double(EEGs.data(ch,:,idxRLD)),3)); RLD = RLD(:)';

                    LS = mean([FLS; RLS], 1, 'omitnan');
                    LD = mean([FLD; RLD], 1, 'omitnan');
                    GS = mean([FLS; FLD], 1, 'omitnan');
                    GD = mean([RLS; RLD], 1, 'omitnan');

                    tt = double(EEGs.times(:))';
                    winIdx = tt >= tStart & tt <= tEnd;

                    if eff == "Loc"
                        stdCurve = LS;
                        devCurve = LD;
                    else
                        stdCurve = GS;
                        devCurve = GD;
                    end

                    if any(winIdx)
                        yStdWin(s) = mean(stdCurve(winIdx), 'omitnan');
                        yDevWin(s) = mean(devCurve(winIdx), 'omitnan');
                    end

                    if hasROIsig
                        [~, ixPeak] = min(abs(tt - tPeakROI));
                        yStdPeak(s) = stdCurve(ixPeak);
                        yDevPeak(s) = devCurve(ixPeak);
                    end
                end

                % ---------------- WINDOW-AVERAGED ROI SUMMARY ----------------
                m = isfinite(yStdWin) & isfinite(yDevWin);
                a = yStdWin(m);
                b = yDevWin(m);
                dff = b - a;   % Deviant - Standard
                n = numel(dff);

                fprintf('\n%s | window-averaged summary\n', roi);
                fprintf('  ROI has significant sample in this window: %s\n', string(hasROIsig));

                if n < 2
                    fprintf('  Not enough valid pairs.\n');
                else
                    meanStd = mean(a);
                    seStd   = std(a,0) / sqrt(n);

                    meanDev = mean(b);
                    seDev   = std(b,0) / sqrt(n);

                    meanDiff = mean(dff);
                    seDiff   = std(dff,0) / sqrt(n);

                    [~, pWin, ciWin, statsWin] = ttest(b, a);
                    dzWin = meanDiff / std(dff,0);

                    fprintf('  Standard mean ± SE = %.4f ± %.4f uV\n', meanStd, seStd);
                    fprintf('  Deviant  mean ± SE = %.4f ± %.4f uV\n', meanDev, seDev);
                    fprintf('  Deviant - Standard mean ± SE = %.4f ± %.4f uV\n', meanDiff, seDiff);
                    fprintf('  Paired t-test on window means: t(%d) = %.4f, p = %.8f, dz = %.4f\n', ...
                        statsWin.df, statsWin.tstat, pWin, dzWin);
                    fprintf('  95%% CI of mean difference = [%.4f, %.4f] uV\n', ciWin(1), ciWin(2));

                    if hasROIsig
                        fprintf('  Minimum corrected pFWE inside ROI/window = %.8f at %.0f ms\n', pminROI, tPeakROI);
                    else
                        fprintf('  No ROI-specific corrected sample inside this scalp-level window.\n');
                    end
                end

                % ---------------- PEAK-SAMPLE ROI SUMMARY ----------------
                if hasROIsig
                    m = isfinite(yStdPeak) & isfinite(yDevPeak);
                    a = yStdPeak(m);
                    b = yDevPeak(m);
                    dff = b - a;
                    n = numel(dff);

                    fprintf('%s | peak-sample summary at %.0f ms\n', roi, tPeakROI);

                    if n < 2
                        fprintf('  Not enough valid peak pairs.\n');
                    else
                        meanStd = mean(a);
                        seStd   = std(a,0) / sqrt(n);

                        meanDev = mean(b);
                        seDev   = std(b,0) / sqrt(n);

                        meanDiff = mean(dff);
                        seDiff   = std(dff,0) / sqrt(n);

                        [~, pPeak, ciPeak, statsPeak] = ttest(b, a);
                        dzPeak = meanDiff / std(dff,0);

                        fprintf('  Standard mean ± SE = %.4f ± %.4f uV\n', meanStd, seStd);
                        fprintf('  Deviant  mean ± SE = %.4f ± %.4f uV\n', meanDev, seDev);
                        fprintf('  Deviant - Standard mean ± SE = %.4f ± %.4f uV\n', meanDiff, seDiff);
                        fprintf('  Paired t-test at peak sample: t(%d) = %.4f, p = %.8f, dz = %.4f\n', ...
                            statsPeak.df, statsPeak.tstat, pPeak, dzPeak);
                        fprintf('  95%% CI of mean difference = [%.4f, %.4f] uV\n', ciPeak(1), ciPeak(2));
                        fprintf('  Corrected pFWE at this ROI/sample = %.8f\n', pminROI);
                    end
                end
            end
        end
    end
end


%% ========================= LOCAL FUNCTIONS =========================

function EEGout = reorder_to_template(EEGin, templateLabels)
curLabels = {EEGin.chanlocs.labels};
[~, idx] = ismember(templateLabels, curLabels);
if any(idx==0)
    missing = templateLabels(idx==0);
    error('Missing channels: %s', strjoin(missing, ', '));
end
EEGout = EEGin;
EEGout.data     = EEGin.data(idx,:,:);
EEGout.chanlocs = EEGin.chanlocs(idx);
EEGout.nbchan   = numel(idx);
EEGout = eeg_checkset(EEGout);
end

function S = to_strlist(v)
if isstring(v)
    S = v;
elseif ischar(v)
    S = string(v);
elseif iscell(v)
    S = strings(0,1);
    for i=1:numel(v)
        if isstring(v{i}) || ischar(v{i})
            S(end+1,1) = string(v{i}); %#ok<AGROW>
        end
    end
else
    S = strings(0,1);
end
S = strtrim(regexprep(S,'\s+',' '));
end

function [devData, stdData] = build_balanced_trials(DATA, devA, devB, stdA, stdB)
devA = devA(:)'; devB = devB(:)'; stdA = stdA(:)'; stdB = stdB(:)';

nDevMin = min(numel(devA), numel(devB));
nStdMin = min(numel(stdA), numel(stdB));

devA = devA(randperm(numel(devA), nDevMin));
devB = devB(randperm(numel(devB), nDevMin));
stdA = stdA(randperm(numel(stdA), nStdMin));
stdB = stdB(randperm(numel(stdB), nStdMin));

devIdx = [devA devB];
stdIdx = [stdA stdB];

nMin = min(numel(devIdx), numel(stdIdx));
devIdx = devIdx(randperm(numel(devIdx), nMin));
stdIdx = stdIdx(randperm(numel(stdIdx), nMin));

devData = double(DATA(:,:,devIdx));
stdData = double(DATA(:,:,stdIdx));
end

function p = unpaired_pmap(devData, stdData, mode)
% devData/stdData: [chan x time x trials]
% mode: "student" (equal-var) or "welch" (unequal-var)
devData = double(devData);
stdData = double(stdData);

n1 = size(devData,3);
n2 = size(stdData,3);

m1 = mean(devData, 3);
m2 = mean(stdData, 3);

v1 = var(devData, 0, 3);
v2 = var(stdData, 0, 3);

if strcmpi(mode,"welch")
    se2 = v1./n1 + v2./n2;
    t   = (m1 - m2) ./ sqrt(se2);
    df  = (se2.^2) ./ ( (v1.^2)./(n1^2*(n1-1)) + (v2.^2)./(n2^2*(n2-1)) );
else
    % Student t-test with pooled variance (paper-safe default)
    sp2 = ((n1-1).*v1 + (n2-1).*v2) ./ (n1+n2-2);
    se2 = sp2 .* (1/n1 + 1/n2);
    t   = (m1 - m2) ./ sqrt(se2);
    df  = (n1+n2-2) * ones(size(t));
end

p = 2 * tcdf(-abs(t), df);
p(~isfinite(p)) = 1;
p = max(min(p,1),0);
end

function p = paired_pmap(devSub, stdSub)
% devSub/stdSub: [sub x chan x time]
devSub = double(devSub);
stdSub = double(stdSub);

D = devSub - stdSub;
n = size(D,1);

m = squeeze(mean(D,1));
s = squeeze(std(D,0,1));

t  = m ./ (s ./ sqrt(n));
df = n - 1;

p = 2 * tcdf(-abs(t), df);
p(~isfinite(p)) = 1;
p = max(min(p,1),0);
end

function pAdj = fdr_bh_map(pCT, q)
% BH-FDR across ALL chan×time entries
p = pCT(:);
[ps, idx] = sort(p,'ascend');
m = numel(ps);
th = (1:m)'/m * q;
k = find(ps <= th, 1, 'last');

pAdj = ones(size(p));
if ~isempty(k)
    adj = ps .* m ./ (1:m)';
    adj = flipud(cummin(flipud(adj)));
    adj = min(adj,1);
    pAdj(idx) = adj;
end
pAdj = reshape(pAdj, size(pCT));
end

function maskTime = triple_criterion_time_mask(pCT, alpha, minElec, minConsec)
% "Time is significant" only if >=minElec electrodes have p<=alpha,
% and this condition lasts >=minConsec consecutive samples.
sigElecCount = sum(pCT <= alpha, 1);
mask0 = sigElecCount >= minElec;
maskTime = false(size(mask0));

d = diff([false mask0 false]);
runStarts = find(d==1);
runEnds   = find(d==-1)-1;
for i = 1:numel(runStarts)
    if (runEnds(i)-runStarts(i)+1) >= minConsec
        maskTime(runStarts(i):runEnds(i)) = true;
    end
end
end

function idxLvl = p_to_binned_index_timecourse(pCT, maskTime, pBins, minElec)
% pCT: [chan x time]
% idxLvl(t) = 0 if not in extent mask; else bin according to pBest(t)=min_elec p

pSorted = sort(pCT, 1, 'ascend');
pBest   = pSorted(minElec, :);
%pBest = min(pCT, [], 1);

T = numel(maskTime);
idxLvl = zeros(1, T);

for t = 1:T
    if ~maskTime(t)
        idxLvl(t) = 0;
    else
        pv = pBest(t);
        lvl = 0;
        for k = 1:numel(pBins)
            if pv <= pBins(k)
                lvl = k;
            end
        end
        idxLvl(t) = lvl; % 0..numel(pBins)
    end
end
end

function plot_erps_2x2(EEGtmpl, tms, STD, DEV, pMap, maskTime, cfg, thisExp, whichEffect, nSub, nStdTrials, nDevTrials, outPng, pMapInt, maskTimeInt)

if nargin < 15
    pMapInt = [];
    maskTimeInt = [];
end

chLabels = {EEGtmpl.chanlocs.labels};

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for ci = 1:4
    chanName = cfg.erpChans{ci};
    ch = find(strcmpi(chLabels, chanName), 1);
    if isempty(ch), ch = 1; end

    nexttile; hold on; grid on;

    mStd = squeeze(mean(STD(:,ch,:),1));
    mDev = squeeze(mean(DEV(:,ch,:),1));

    semStd = squeeze(std(STD(:,ch,:),0,1)) ./ sqrt(size(STD,1));
    semDev = squeeze(std(DEV(:,ch,:),0,1)) ./ sqrt(size(DEV,1));

    draw_sem_fill(tms, mStd, semStd, [0 0.45 0.85], 0.15);
    draw_sem_fill(tms, mDev, semDev, [0.90 0.45 0.10], 0.15);

    plot(tms, mStd, 'LineWidth', 2, 'Color', [0 0.45 0.85]);
    plot(tms, mDev, 'LineWidth', 2, 'Color', [0.90 0.45 0.10]);

    xline(0,'k-','LineWidth',1);
    yline(0,'k-','LineWidth',1);

    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    ylim(cfg.ylim);
    xlim([tms(1) tms(end)]);

    title(sprintf('ERP at %s', upper(chanName)), 'FontWeight','bold');

    if strcmpi(whichEffect,'LOCAL')
        legStd = 'Local standard';
        legDev = 'Local deviant';
    elseif strcmpi(whichEffect,'GLOBAL')
        legStd = 'Global standard';
        legDev = 'Global deviant';
    else
        legStd = 'Standard';
        legDev = 'Deviant';
    end

    legend({legStd, legDev}, 'Location','southwest');

    tmsR  = tms(:)';
    maskR = logical(maskTime(:)');

    pch = pMap(ch,:);
    pch = pch(:)';

    showMask = maskR & (tmsR >= cfg.showFromMs) & (pch <= cfg.alpha_grp);
    showMask = bridge_small_gaps(showMask, cfg.barMaxGap);
    showMask = keep_min_run_length(showMask, cfg.barMinConsec);

    % ================= BLUE main-effect bar =================
    if any(showMask)
        yBar = cfg.ylim(2) - cfg.barYoffset_uV;

        d = diff([false showMask false]);
        runStarts = find(d == 1);
        runEnds   = find(d == -1) - 1;

        for rr = 1:numel(runStarts)
            idx = runStarts(rr):runEnds(rr);
            plot(tmsR(idx), yBar*ones(size(idx)), '-', ...
                'LineWidth', 6, 'Color', [0 0 1], 'HandleVisibility','off');
        end

        if isfield(cfg,'pText') && cfg.pText && ~isempty(runStarts)
            firstIdx = runStarts(1);

            useMCC = isfield(cfg,'mcc') && isfield(cfg.mcc,'method') && ...
                ~strcmpi(string(cfg.mcc.method),'none') && should_apply_mcc(cfg,'group');

            if useMCC
                lab = sprintf('pFWE<%.2f', cfg.alpha_grp);
            else
                lab = sprintf('p<%.2f', cfg.alpha_grp);
            end

            text(gca, tmsR(firstIdx), yBar+0.8, lab, ...
                'Color', [0 0 1], 'FontSize', 9, 'FontWeight','bold');
        end
    end

    % ================= GREY interaction bar =================
    if ~isempty(pMapInt) && ~isempty(maskTimeInt)
        pchInt   = pMapInt(ch,:);
        pchInt   = pchInt(:)';
        maskIntR = logical(maskTimeInt(:)');

        showMaskInt = maskIntR & (tmsR >= cfg.showFromMs) & (pchInt <= cfg.alpha_grp);
        showMaskInt = bridge_small_gaps(showMaskInt, cfg.barMaxGap);
        showMaskInt = keep_min_run_length(showMaskInt, cfg.barMinConsec);

        if any(showMaskInt)
            yBarInt = cfg.ylim(2) - cfg.barYoffset_uV - 0.55;

            dInt = diff([false showMaskInt false]);
            runStartsInt = find(dInt == 1);
            runEndsInt   = find(dInt == -1) - 1;

            for rr = 1:numel(runStartsInt)
                idx = runStartsInt(rr):runEndsInt(rr);
                plot(tmsR(idx), yBarInt*ones(size(idx)), '-', ...
                    'LineWidth', 5, 'Color', [0.45 0.45 0.45], 'HandleVisibility','off');
            end
        end
    end
end

if strcmpi(thisExp,'E1')
    sgtitle_experiment = 'Counting';
else
    sgtitle_experiment = 'Visual interference';
end

sgtitle(sprintf('Event Related Potentials | %s | %s effects', ...
    sgtitle_experiment, [upper(whichEffect(1)), lower(whichEffect(2:end))]), ...
    'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function plot_topo_snaps_erp(EEGtmpl, tms, diffCT, cfg, bigTitle, outPng)

snap   = cfg.snapTimes_ms(:)';
chlocs = EEGtmpl.chanlocs;
nSnap  = numel(snap);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.38]);

leftMargin   = 0.012;
rightMargin  = 0.012;
topMargin    = 0.16;
bottomMargin = 0.20;
gapX         = 0;
usableW = 1 - leftMargin - rightMargin - (nSnap-1)*gapX;
axW = usableW / nSnap;
axH = 0.58;
y0  = bottomMargin + 0.08;

ax = gobjects(nSnap,1);

for i = 1:nSnap
    x0 = leftMargin + (i-1)*(axW + gapX);

    ax(i) = axes('Parent',fig,'Units','normalized', ...
        'Position',[x0 y0 axW axH]);

    [~, ti] = min(abs(tms - snap(i)));

    topoplot(double(diffCT(:,ti)), chlocs, ...
        'electrodes','on', ...
        'numcontour', cfg.topoContours, ...
        'gridscale',  cfg.topoGrid, ...
        'shading',    cfg.topoShading, ...
        'whitebk',    'off');

    caxis(ax(i), cfg.caxis_erpTopo);

    title(ax(i), sprintf('%d ms', snap(i)), ...
        'FontSize', 12, 'FontWeight','bold');
end

sgtitle(fig, bigTitle, 'FontWeight','bold');

cb = colorbar('southoutside');
cb.Label.String = '\muV (Dev-Std)';
cb.Label.FontSize = 12;
cb.Position = [0.37 0.11 0.25 0.045];

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function plot_topo_snaps_stats(EEGtmpl, tms, pCT, maskTime, cfg, bigTitle, outPng)
% plot_topo_snaps_stats
% One-row compact stats topoplots using the SAME layout logic as
% eeg_plot_topo_snaps_condition, but plotting p-values and red markers.

snap   = cfg.snapTimes_ms(:)';
chlocs = EEGtmpl.chanlocs;
nSnap  = numel(snap);

% --- RP-style blue gradient ---
a = [247,251,255]/255;
b = [222,235,247]/255;
c = [198,219,239]/255;
d = [158,202,225]/255;
e = [107,174,214]/255;
f = [66,146,198]/255;
g = [33,113,181]/255;
customColors   = [g; f; e; d; c; b; a];
smoothGradient = interp1(linspace(0,1,size(customColors,1)), ...
    customColors, linspace(0,1,256), 'pchip');

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.38]);

% ---------- SAME manual compact layout as your working condition topos ----------
leftMargin   = 0.012;
rightMargin  = 0.012;
topMargin    = 0.16;
bottomMargin = 0.20;   % leave room for colorbar
gapX         = 0;  % smaller = tighter
usableW = 1 - leftMargin - rightMargin - (nSnap-1)*gapX;
axW = usableW / nSnap;
axH = 0.58;
y0  = bottomMargin + 0.08;

ax = gobjects(nSnap,1);

for i = 1:nSnap
    x0 = leftMargin + (i-1)*(axW + gapX);

    ax(i) = axes('Parent',fig,'Units','normalized', ...
        'Position',[x0 y0 axW axH]);

    [~, ti] = min(abs(tms - snap(i)));

    sigCh = [];
    if maskTime(ti)
        sigCh = find(pCT(:,ti) <= cfg.alpha_grp);
    end

    if isempty(sigCh)
        topoplot(double(pCT(:,ti)), chlocs, ...
            'electrodes','off', ...
            'maplimits',[0 1], ...
            'numcontour', cfg.topoContours, ...
            'gridscale',  cfg.topoGrid, ...
            'shading',    cfg.topoShading_stats, ...
            'whitebk',    'off');
    else
        topoplot(double(pCT(:,ti)), chlocs, ...
            'electrodes','off', ...
            'maplimits',[0 1], ...
            'numcontour', cfg.topoContours, ...
            'gridscale',  cfg.topoGrid, ...
            'shading',    cfg.topoShading_stats, ...
            'whitebk',    'off', ...
            'emarker2',{sigCh,'o','r',5,0.8});   % smaller red dots
    end

    title(ax(i), sprintf('%d ms', snap(i)), ...
        'FontSize', 12, 'FontWeight','bold');
end

colormap(fig, smoothGradient);

sgtitle(fig, sprintf('%s', ...
    bigTitle), 'FontWeight','bold');

% ---------- SAME style of colorbar placement ----------
axCB = axes('Parent',fig,'Units','normalized', ...
    'Position',[0 0 1 1], ...
    'Visible','off', ...
    'HandleVisibility','off');
colormap(axCB, smoothGradient);
set(axCB,'CLim',[0 1]);

cb = colorbar(axCB,'southoutside');
cb.Label.String = 'p-Values';
cb.Label.FontSize = 12;
cb.Position = [0.37 0.11 0.25 0.045];

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end


function plot_subjectwise_matrix_post0(SSI, cfg, outPng)
% Subject-wise: subjects × time, value is binned by pBest(t) with extent mask applied
tms = SSI.E1.times;
keep = tms >= 0;
tms = tms(keep);

E1L = SSI.E1.local(:,keep);
E1G = SSI.E1.global(:,keep);
E2L = SSI.E2.local(:,keep);
E2G = SSI.E2.global(:,keep);

% Optional display hiding (rarely needed)
if isfinite(cfg.displayHideBinsAboveP)
    % hide everything weaker than that p (i.e., bins where threshold > desired)
    pBins = cfg.pBins_matrix;
    hideLvl = find(pBins <= cfg.displayHideBinsAboveP, 1, 'last');
    if isempty(hideLvl), hideLvl = 0; end
    E1L(E1L < hideLvl) = 0; E1G(E1G < hideLvl) = 0;
    E2L(E2L < hideLvl) = 0; E2G(E2G < hideLvl) = 0;
end

% Colormap + ticks
cmap = pbin_cmap(cfg);
nLvl = numel(cfg.pBins_matrix);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

% ===== FIX: leave space at the right for TWO colorbars =====
t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
t.OuterPosition = [0.03 0.04 0.88 0.92]; % [left bottom width height] leaving right margin

nexttile;
imagesc(tms, 1:size(E1L,1), E1L);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants'); % << CHANGED
title('Counting | Local effects');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

nexttile;
imagesc(tms, 1:size(E1G,1), E1G);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants'); % << CHANGED
title('Counting | Global effects');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

nexttile;
imagesc(tms, 1:size(E2L,1), E2L);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants'); % << CHANGED
title('Visual interference | Local effects');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

nexttile;
imagesc(tms, 1:size(E2G,1), E2G);
xlabel('Time after 5th sound onset (ms)'); ylabel('Participants'); % << CHANGED
title('Visual interference | Global effects');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

% ===== FIX: TWO identical colorbars (top row + bottom row) =====
% Pass label preferences via fig.UserData (robust without changing signatures everywhere)
fig.UserData.noneLabel = char(getfield_default(cfg,'noneLabel', '>0.05')); %#ok<GFLD>
labelStr = char(getfield_default(cfg,'binLabel', 'p-values'));            %#ok<GFLD>

add_two_row_colorbars_right(fig, cfg.pBins_matrix, cmap, labelStr, [-0.5 nLvl+0.5]);

% sgtitle(sprintf('Bekinschtein-like SUBJECT matrix (subjects×time) | minElec=%d, minConsec=%d | (t>=0 only)', ...
%     cfg.minElec_ind, cfg.minConsec_ind), 'FontWeight','bold');

sgtitle('Local and global effects in individual subjects | Event Related Potentials', 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function plot_mean_subjectwise_index_post0(SSI, cfg, outPng)
tms = SSI.E1.times;
keep = tms >= 0;
tms = tms(keep);

mE1L = mean(SSI.E1.local(:,keep),1);
mE1G = mean(SSI.E1.global(:,keep),1);
mE2L = mean(SSI.E2.local(:,keep),1);
mE2G = mean(SSI.E2.global(:,keep),1);

% Auto y-limit same for both panels
yMax = max([mE1L(:); mE1G(:); mE2L(:); mE2G(:)]);
yMax = yMax * 1.05;
if yMax < 1, yMax = 1; end

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

% =================== Counting ===================
nexttile; hold on; grid on;
plot(tms, mE1L, 'LineWidth', 2, 'Color',[0 0.45 0.85]);
plot(tms, mE1G, 'LineWidth', 2, 'Color',[0.90 0.45 0.10]);
xlabel('Time after 5th sound onset (ms)');
ylabel('Mean binned index');
title('Counting');
legend({'Local effect','Global effect'}, 'Location','northeast');
ylim([0 yMax]);

% =================== Visual interference ===================
nexttile; hold on; grid on;
plot(tms, mE2L, 'LineWidth', 2,'Color',[0 0.45 0.85]);
plot(tms, mE2G, 'LineWidth', 2,'Color',[0.90 0.45 0.10]);
xlabel('Time after 5th sound onset (ms)');
ylabel('Mean binned index');
title('Visual interference');
legend({'Local effect','Global effect'}, 'Location','northeast');
ylim([0 yMax]);

%sgtitle('Mean subject-wise binned index | post-0 only', 'FontWeight','bold');
sgtitle('Mean individual statistical significance index | Event Related Potentials', 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function plot_group_pval_matrix_2x2(GROUPPMAP, cfg, outPng)
% Electrode-wise (group): electrodes × time, group paired t-test p-values,
% with extent rules applied by maskTime (>=minElec for >=minConsec).

tms = GROUPPMAP.E1.times;
keep = tms >= 0;
tms = tms(keep);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

% ===== FIX: leave space at the right for TWO colorbars =====
t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
t.OuterPosition = [0.03 0.04 0.88 0.92]; % leaving right margin

pb = cfg.pBins_matrix;     % bins used for display
nLvl = numel(pb);
cmap = pbin_cmap(cfg);     % consistent across BOTH figures

    function M = makeMat(pCT, maskTime)
        pCT = pCT(:,keep);
        maskTime = maskTime(keep);

        M = zeros(size(pCT)); % 0 = none
        for tt = 1:numel(maskTime)
            if ~maskTime(tt), continue; end
            pv = pCT(:,tt);
            lvl = zeros(size(pv));
            for k=1:numel(pb)
                lvl(pv <= pb(k)) = k; % higher k = stricter p
            end
            M(:,tt) = lvl;
        end
    end

M_E1L = makeMat(GROUPPMAP.E1.pLoc, GROUPPMAP.E1.mLoc);
M_E1G = makeMat(GROUPPMAP.E1.pGlo, GROUPPMAP.E1.mGlo);
M_E2L = makeMat(GROUPPMAP.E2.pLoc, GROUPPMAP.E2.mLoc);
M_E2G = makeMat(GROUPPMAP.E2.pGlo, GROUPPMAP.E2.mGlo);

nexttile;
imagesc(tms, 1:size(M_E1L,1), M_E1L);
xlabel('Time after 5th sound onset (ms)'); ylabel('Electrodes'); % << CHANGED
title('Counting | Local effects');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

nexttile;
imagesc(tms, 1:size(M_E1G,1), M_E1G);
xlabel('Time after 5th sound onset (ms)'); ylabel('Electrodes'); % << CHANGED
title('Counting | Global effects');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

nexttile;
imagesc(tms, 1:size(M_E2L,1), M_E2L);
xlabel('Time after 5th sound onset (ms)'); ylabel('Electrodes'); % << CHANGED
title('Visual interference | Local');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

nexttile;
imagesc(tms, 1:size(M_E2G,1), M_E2G);
xlabel('Time after 5th sound onset (ms)'); ylabel('Electrodes'); % << CHANGED
title('Visual interference | Global');
colormap(gca, cmap); caxis([-0.5 nLvl+0.5]);
% --- FIX: ensure x-axis shows 0 and 700 (epoch limit) ---
xlim([cfg.tlim_stats(1) cfg.tlim_stats(2)]);
xt = get(gca,'XTick');
set(gca,'XTick', unique([cfg.tlim_stats(1) xt cfg.tlim_stats(2)]));

% ===== FIX: TWO identical colorbars (top row + bottom row) =====
fig.UserData.noneLabel = char(getfield_default(cfg,'noneLabel', '>0.05')); %#ok<GFLD>
labelStr = char(getfield_default(cfg,'binLabel', 'p-values'));            %#ok<GFLD>

add_two_row_colorbars_right(fig, cfg.pBins_matrix, cmap, labelStr, [-0.5 nLvl+0.5]);

% If MC was enabled, mention it in the title (optional)
if isfield(cfg,'mc') && isfield(cfg.mc,'doMC') && cfg.mc.doMC
    extra = sprintf(' | MC-validated (alpha=%.2f, nPerm=%d)', cfg.mc.alpha, cfg.mc.nPerm);
else
    extra = '';
end

sgtitle(['Bekinschtein-like GROUP significance matrix (electrodes×time; extent rules; t>=0 ms)' extra], 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function cmap = pbin_cmap(cfg)
% returns (nBins+1) x 3 colormap for levels: 0..nBins

nBins = numel(cfg.pBins_matrix);
needN = nBins + 1;

% Preferred (new): cfg.binColors
if isfield(cfg,'binColors') && ~isempty(cfg.binColors)
    C = cfg.binColors;
    if size(C,2) ~= 3 || size(C,1) ~= needN
        error('cfg.binColors must be [%d x 3] (nBins+1).', needN);
    end
    cmap = C;
    return;
end

% Fallback (old): cfg.pBinColors (kept for backward compatibility)
cols = nan(needN,3);
cols(1,:) = cfg.pBinColors.none;

% Expect fields c1..cN
for k = 1:nBins
    fn = sprintf('c%d', k);
    if ~isfield(cfg.pBinColors, fn)
        error('Missing cfg.pBinColors.%s for %d bins.', fn, nBins);
    end
    cols(k+1,:) = cfg.pBinColors.(fn);
end

cmap = cols;
end

function [ticks, labels] = pbin_ticks_labels(pBins, noneLabel)
% 0..nBins with a custom label at 0 (e.g. ">0.05")
n = numel(pBins);
ticks = 0:n;
labels = cell(1,n+1);

if nargin < 2 || isempty(noneLabel)
    noneLabel = '>0.05';
end
labels{1} = char(noneLabel);

for k = 1:n
    labels{k+1} = sprintf('<=%.4g', pBins(k));
end
end

function add_two_row_colorbars_right(fig, pBins, cmap, labelStr, clim)
% Adds TWO identical discrete colorbars on the right side of the figure:
% one aligned with the TOP row, one aligned with the BOTTOM row.
%
% Works across MATLAB versions by attaching each colorbar to a dummy axes.

noneLabel = '>0.05';
if isfield(fig.UserData,'noneLabel')
    noneLabel = fig.UserData.noneLabel;
end
[ticks, labels] = pbin_ticks_labels(pBins, noneLabel);

% --- Top dummy axes + colorbar ---
axTop = axes('Parent', fig, 'Units','normalized', 'Position',[0.915 0.56 0.02 0.36], ...
    'Visible','off', 'HandleVisibility','off');
colormap(axTop, cmap);
set(axTop, 'CLim', clim);

cb1 = colorbar(axTop, 'Units','normalized');
cb1.Position = [0.94 0.56 0.015 0.36]; % [x y w h]
cb1.Ticks = ticks;
cb1.TickLabels = labels;
cb1.Label.String = labelStr;

% --- Bottom dummy axes + colorbar ---
axBot = axes('Parent', fig, 'Units','normalized', 'Position',[0.915 0.10 0.02 0.36], ...
    'Visible','off', 'HandleVisibility','off');
colormap(axBot, cmap);
set(axBot, 'CLim', clim);

cb2 = colorbar(axBot, 'Units','normalized');
cb2.Position = [0.94 0.10 0.015 0.36]; % [x y w h]
cb2.Ticks = ticks;
cb2.TickLabels = labels;
cb2.Label.String = labelStr;
end

function mask = keep_min_run_length(mask, minLen)
mask = mask(:)';
d = diff([false mask false]);
rs = find(d==1);
re = find(d==-1)-1;
out = false(size(mask));
for i=1:numel(rs)
    if (re(i)-rs(i)+1) >= minLen
        out(rs(i):re(i)) = true;
    end
end
mask = out;
end

function mask = bridge_small_gaps(mask, maxGap)
mask = mask(:)';
if maxGap <= 0, return; end
x = mask;

d = diff([false x false]);
rs = find(d==1);
re = find(d==-1)-1;

if isempty(rs), return; end

for i=1:numel(rs)-1
    gapStart = re(i)+1;
    gapEnd   = rs(i+1)-1;
    gapLen   = gapEnd-gapStart+1;
    if gapLen > 0 && gapLen <= maxGap
        x(gapStart:gapEnd) = true;
    end
end
mask = x;
end

function expand_tiled_axes(ax, left, bottom, width, height)
ax = ax(isgraphics(ax));
if isempty(ax), return; end

pos = cell2mat(get(ax,'Position'));
minX = min(pos(:,1)); minY = min(pos(:,2));
maxX = max(pos(:,1)+pos(:,3));
maxY = max(pos(:,2)+pos(:,4));

curW = maxX - minX;
curH = maxY - minY;

for i=1:size(pos,1)
    x = (pos(i,1)-minX)/curW;
    y = (pos(i,2)-minY)/curH;
    w = pos(i,3)/curW;
    h = pos(i,4)/curH;

    newPos = [left + x*width, bottom + y*height, w*width, h*height];
    set(ax(i),'Position',newPos);
end
end

% ===================== NEW: MONTE-CARLO HELPERS (GROUP, paired) =====================

function [pMC, nHit] = mc_validate_triplecriterion_group(devSub, stdSub, alpha, minElec, minConsec, nPerm, seed, verbose, progressStep)
% mc_validate_triplecriterion_group
% Paper-like Monte-Carlo validation of the GROUP triple-criterion.
%
% Inputs:
%   devSub/stdSub : [sub x chan x time] subject-level ERPs (already cropped to stats window)
%   alpha/minElec/minConsec : same as the triple-criterion parameters
%   nPerm, seed : Monte-Carlo settings
%
% Output:
%   pMC  : (1 + nHit) / (1 + nPerm)
%   nHit : # permutations where the triple-criterion occurs anywhere in time
%
% Null model (paired):
%   randomly swap dev/std within each subject (equivalent to sign-flip in D=dev-std).
%
% Criterion counted:
%   whether there exists ANY run of >=minConsec time samples where >=minElec electrodes have p<=alpha.

rng(seed);

% Observed criterion (for debugging; not strictly needed for the paper-like pMC)
% pObs = paired_pmap(devSub, stdSub);
% obsHasCrit = triple_criterion_exists_anywhere(pObs, alpha, minElec, minConsec);

nHit = 0;

nSub = size(devSub,1);
if nSub < 2
    pMC = 1;
    nHit = nPerm;
    return;
end

for ip = 1:nPerm
    % Random swap mask per subject
    swap = rand(nSub,1) > 0.5;

    devP = devSub;
    stdP = stdSub;

    % Swap dev/std for chosen subjects
    if any(swap)
        devP(swap,:,:) = stdSub(swap,:,:);
        stdP(swap,:,:) = devSub(swap,:,:);
    end

    % Paired t-test p-map under permutation
    pPerm = paired_pmap(devP, stdP);

    % Does the extent criterion happen anywhere?
    hasCrit = triple_criterion_exists_anywhere(pPerm, alpha, minElec, minConsec);

    if hasCrit
        nHit = nHit + 1;
    end

    if verbose && mod(ip, progressStep) == 0
        fprintf('      MC perm %d/%d | hits=%d\n', ip, nPerm, nHit);
    end
end

% Monte-Carlo p-value (add-one smoothing)
pMC = (1 + nHit) / (1 + nPerm);
end

function tf = triple_criterion_exists_anywhere(pCT, alpha, minElec, minConsec)
% Returns true if there is ANY time-run of length >=minConsec where
% >=minElec electrodes have p<=alpha.
sigElecCount = sum(pCT <= alpha, 1);
mask0 = sigElecCount >= minElec;

d = diff([false mask0 false]);
rs = find(d==1);
re = find(d==-1)-1;

tf = false;
for i = 1:numel(rs)
    if (re(i)-rs(i)+1) >= minConsec
        tf = true;
        return;
    end
end
end

% ---- helper: safe getfield default ----
function v = getfield_default(S, field, defaultVal)
if isstruct(S) && isfield(S, field)
    v = S.(field);
else
    v = defaultVal;
end
end

function [pAdj, info] = apply_mcc_pmap(pCT, cfg, tag)
% apply_mcc_pmap
% Generic MCC for a p-map (chan x time).
% Methods:
%   none, fdr_bh, holm, bonferroni
% NOTE: fwe_permutation is NOT supported here (needs test-stat + permutation design).

info = struct();
if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'method') || isempty(cfg.mcc.method)
    info.method = 'none';
    pAdj = pCT;
    return;
end

method = lower(string(cfg.mcc.method));
info.method = char(method);

switch method
    case "none"
        pAdj = pCT;

    case "fdr_bh"
        q = getfield_default(cfg.mcc,'q', 0.05);
        pAdj = fdr_bh_map(pCT, q);
        info.q = q;

    case "holm"
        pAdj = holm_bonferroni_map(pCT);

    case "bonferroni"
        pAdj = bonferroni_map(pCT);

    case "fwe_permutation"
        % Not valid without a permutation distribution for a specific test.
        warning('[%s] MCC method=fwe_permutation not supported for generic p-maps. Skipping (no correction).', tag);
        pAdj = pCT;

    otherwise
        warning('[%s] Unknown MCC method: %s. Using none.', tag, method);
        pAdj = pCT;
        info.method = 'none';
end
end

function [pAdj, info] = apply_mcc_group(devSub, stdSub, pAlready, cfg, tag)
% apply_mcc_group
% MCC for GROUP paired maps.
% If cfg.mcc.method == 'fwe_permutation', compute max-|t| permutation FWE p-map.
% Otherwise, apply p-map corrections (FDR/Holm/Bonferroni/none).

info = struct();
if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'method') || isempty(cfg.mcc.method)
    info.method = 'none';
    pAdj = pAlready;
    return;
end

method = lower(string(cfg.mcc.method));
info.method = char(method);

switch method
    case "fwe_permutation"
        % Paired sign-flip max-|t| FWE
        fwe = cfg.mcc.fwe;
        nPerm = getfield_default(fwe,'nPerm',2000);
        seed  = getfield_default(fwe,'seed',777);
        verbose = getfield_default(fwe,'verbose',true);
        step  = getfield_default(fwe,'progressStep',200);
        [pAdj, info] = fwe_maxT_paired(devSub, stdSub, nPerm, seed, verbose, step);
        info.method = sprintf('fwe_permutation(maxT, nPerm=%d)', nPerm);

    otherwise
        % Apply generic p-map MCC
        [pAdj, info2] = apply_mcc_pmap(pAlready, cfg, tag);
        info = info2;
end
end

function pAdj = bonferroni_map(pCT)
% Bonferroni correction across ALL chan×time entries.
p = pCT(:);
m = numel(p);
pAdj = min(p * m, 1);
pAdj = reshape(pAdj, size(pCT));
end

function pAdj = holm_bonferroni_map(pCT)
% Holm-Bonferroni adjusted p-values (step-down) across ALL chan×time entries.
p = pCT(:);
m = numel(p);
[ps, idx] = sort(p,'ascend');

adj = zeros(size(ps));
for i = 1:m
    adj(i) = (m - i + 1) * ps(i);
end
% Step-down monotonicity: cumulative max
adj = cummax(adj);
adj = min(adj, 1);

pAdj = ones(size(p));
pAdj(idx) = adj;
pAdj = reshape(pAdj, size(pCT));
end

function [pFWE, info] = fwe_maxT_paired(devSub, stdSub, nPerm, seed, verbose, progressStep)
% fwe_maxT_paired
% Permutation FWE correction using max-|t| for paired data via sign flips.
% Inputs:
%   devSub/stdSub: [sub x chan x time]
% Outputs:
%   pFWE: [chan x time] family-wise-error corrected p-values
%   info: struct with details (maxT distribution, etc.)

rng(seed);

% Observed t-map
[tObs, df] = paired_tmap(devSub, stdSub);
absTobs = abs(tObs(:));

% Permutation distribution of max|t|
maxT = zeros(nPerm,1);
nSub = size(devSub,1);

D = double(devSub) - double(stdSub); % [sub x chan x time]

for ip = 1:nPerm
    flip = (rand(nSub,1) > 0.5);
    sgn = ones(nSub,1);
    sgn(flip) = -1;

    Dp = D .* reshape(sgn, [nSub 1 1]);  % sign flip per subject

    % compute t-map for Dp against zero
    m = squeeze(mean(Dp,1));
    s = squeeze(std(Dp,0,1));
    t = m ./ (s ./ sqrt(nSub));
    t(~isfinite(t)) = 0;

    maxT(ip) = max(abs(t(:)));

    if verbose && mod(ip, progressStep) == 0
        fprintf('      FWE maxT perm %d/%d\n', ip, nPerm);
    end
end

% FWE p-values: pFWE(i) = P(max|t_perm| >= |t_obs(i)|)
% Add-one smoothing
maxTsorted = sort(maxT);
pFWE_vec = zeros(size(absTobs));
for i = 1:numel(absTobs)
    % count how many perm maxima >= observed
    nGE = sum(maxTsorted >= absTobs(i));
    pFWE_vec(i) = (1 + nGE) / (1 + nPerm);
end

pFWE = reshape(pFWE_vec, size(tObs));

info = struct();
info.df = df;
info.nPerm = nPerm;
info.seed = seed;
end

function [tMap, df] = paired_tmap(devSub, stdSub)
% paired_tmap
% devSub/stdSub: [sub x chan x time]
D = double(devSub) - double(stdSub);
n = size(D,1);
m = squeeze(mean(D,1));
s = squeeze(std(D,0,1));
tMap = m ./ (s ./ sqrt(n));
tMap(~isfinite(tMap)) = 0;
df = n - 1;
end


function print_run_summary(GROUPPMAP, SSI, cfg)
% print_run_summary
% Prints a human-readable run summary for each experiment (E1/E2):
%   - MCC method (none/fdr_bh/holm/bonferroni/fwe_permutation)
%   - Extent-rule parameters (alpha, minElec, minConsec)
%   - Whether any time survived the extent rule pre-MC and post-MC
%   - If MC enabled: pMC + hit counts, and whether masks were zeroed
%   - Nsub and time window info

fprintf('\n============================================================\n');
fprintf('RUN SUMMARY | Bekinschtein-like EEG stats\n');
fprintf('============================================================\n');

% Time window
if isfield(cfg,'tlim_stats')
    fprintf('Stats time window: [%g %g] ms\n', cfg.tlim_stats(1), cfg.tlim_stats(2));
end

% Extent-rule parameters
fprintf('Extent rule (GROUP):  alpha=%g | minElec=%d | minConsec=%d samples\n', ...
    cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp);
fprintf('Extent rule (INDIV):  alpha=%g | minElec=%d | minConsec=%d samples\n', ...
    cfg.alpha_ind, cfg.minElec_ind, cfg.minConsec_ind);

% Sampling-rate implication (optional, purely informative)
if isfield(cfg,'minConsec_grp') && isfield(cfg,'minConsec_ind')
    fprintf('NOTE: minConsec is in SAMPLES (depends on sampling rate).\n');
end

% MCC method
if isfield(cfg,'mcc') && isfield(cfg.mcc,'method')
    m = char(lower(string(cfg.mcc.method)));
else
    m = 'none';
end
applyTo = 'both';
if isfield(cfg,'mcc') && isfield(cfg.mcc,'applyTo') && ~isempty(cfg.mcc.applyTo)
    applyTo = char(lower(string(cfg.mcc.applyTo)));
end
fprintf('Multiple-comparison correction (MCC): %s | applyTo=%s\n', m, applyTo);

if strcmpi(m,'fdr_bh') && isfield(cfg.mcc,'q')
    fprintf('  FDR q = %.3f\n', cfg.mcc.q);
end
if strcmpi(m,'fwe_permutation') && isfield(cfg.mcc,'fwe')
    fprintf('  FWE perm nPerm=%d | seed=%d\n', ...
        getfield_default(cfg.mcc.fwe,'nPerm',NaN), getfield_default(cfg.mcc.fwe,'seed',NaN));
end

% Monte-Carlo extent validation
doMC = isfield(cfg,'mc') && isfield(cfg.mc,'doMC') && cfg.mc.doMC;
fprintf('Monte-Carlo extent validation (GROUP): %s\n', ternary(doMC,'ON','OFF'));
if doMC
    fprintf('  MC alpha=%.3f | nPerm=%d | zeroMaskIfNotValidated=%s\n', ...
        cfg.mc.alpha, cfg.mc.nPerm, ternary(getfield_default(cfg.mc,'zeroMaskIfNotValidated',true),'true','false'));
end

fprintf('------------------------------------------------------------\n');

exps = intersect(fieldnames(GROUPPMAP), {'E1','E2'}, 'stable');
if isempty(exps)
    fprintf('No GROUPPMAP results found. (Maybe experiments skipped?)\n');
    fprintf('============================================================\n\n');
    return;
end

for i = 1:numel(exps)
    E = exps{i};
    fprintf('%s\n', E);
    mm = getfield_default(GROUPPMAP.(E),'mcc_method','(unknown)');
    aa = getfield_default(GROUPPMAP.(E),'mcc_applyTo','(unknown)');
    fprintf('  MCC used: %s | applyTo=%s', mm, aa);

    if isfield(GROUPPMAP.(E),'mcc_q') && isfinite(GROUPPMAP.(E).mcc_q)
        fprintf(' | q=%.3f', GROUPPMAP.(E).mcc_q);
    end
    if isfield(GROUPPMAP.(E),'mcc_fwe_nPerm') && isfinite(GROUPPMAP.(E).mcc_fwe_nPerm)
        fprintf(' | fwe_nPerm=%d', GROUPPMAP.(E).mcc_fwe_nPerm);
    end
    fprintf('\n');

    % Nsub
    nSub = NaN;
    if isfield(SSI,E) && isfield(SSI.(E),'subIDs')
        nSub = numel(SSI.(E).subIDs);
    end
    if isfinite(nSub)
        fprintf('  Nsub: %d\n', nSub);
    end

    % Raw vs post-MC survival
    hasL_raw  = getfield_default(GROUPPMAP.(E),'obsHasLocal_raw',  any(getfield_default(GROUPPMAP.(E),'mLoc_raw',false)));
    hasG_raw  = getfield_default(GROUPPMAP.(E),'obsHasGlobal_raw', any(getfield_default(GROUPPMAP.(E),'mGlo_raw',false)));
    hasL_post = getfield_default(GROUPPMAP.(E),'obsHasLocal_postMC',  any(getfield_default(GROUPPMAP.(E),'mLoc',false)));
    hasG_post = getfield_default(GROUPPMAP.(E),'obsHasGlobal_postMC', any(getfield_default(GROUPPMAP.(E),'mGlo',false)));

    fprintf('  Extent criterion (pre-MC):  Local=%s | Global=%s\n', ...
        ternary(hasL_raw,'YES','NO'), ternary(hasG_raw,'YES','NO'));

    if doMC
        pMC_L = getfield_default(GROUPPMAP.(E),'pMC_local', NaN);
        pMC_G = getfield_default(GROUPPMAP.(E),'pMC_global', NaN);
        nHitL = getfield_default(GROUPPMAP.(E),'nHit_local', NaN);
        nHitG = getfield_default(GROUPPMAP.(E),'nHit_global', NaN);

        fprintf('  MC p-values: Local=%.4f (hits=%s) | Global=%.4f (hits=%s)\n', ...
            pMC_L, num2str(nHitL), pMC_G, num2str(nHitG));

        fprintf('  Extent criterion (post-MC): Local=%s | Global=%s\n', ...
            ternary(hasL_post,'YES','NO'), ternary(hasG_post,'YES','NO'));

        if getfield_default(cfg.mc,'zeroMaskIfNotValidated',true)
            if hasL_raw && ~hasL_post
                fprintf('    -> Local was ZEROED by MC gating (not validated).\n');
            end
            if hasG_raw && ~hasG_post
                fprintf('    -> Global was ZEROED by MC gating (not validated).\n');
            end
        end
    else
        fprintf('  Extent criterion (post-MC): Local=%s | Global=%s\n', ...
            ternary(hasL_post,'YES','NO'), ternary(hasG_post,'YES','NO'));
    end

    fprintf('------------------------------------------------------------\n');
end

fprintf('============================================================\n\n');

% Optional: estimate sample period if times are available in GROUPPMAP
if isfield(GROUPPMAP.(exps{1}),'times')
    t = GROUPPMAP.(exps{1}).times;
    if numel(t) > 2
        dt = median(diff(t)); % ms
        fprintf('Approx sample period from EEG.times: %.3f ms (fs ~ %.1f Hz)\n', dt, 1000/dt);
        fprintf('GROUP minConsec=%d ~ %.1f ms | IND minConsec=%d ~ %.1f ms\n', ...
            cfg.minConsec_grp, cfg.minConsec_grp*dt, cfg.minConsec_ind, cfg.minConsec_ind*dt);
    end
end

end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


function [pAdj, info] = apply_mcc_individual(devData, stdData, pAlready, cfg, tag, ttestMode)
% apply_mcc_individual
% MCC for INDIVIDUAL unpaired trial-level maps.
% If cfg.mcc.method == 'fwe_permutation', compute max-|t| permutation FWE p-map
% by shuffling condition labels across trials (exchangeability under H0).
%
% devData/stdData: [chan x time x trials]
% pAlready:        [chan x time] raw p-map from unpaired_pmap
%
% NOTE: This is separate from apply_mcc_pmap because permutation FWE needs the DATA.

info = struct();
if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'method') || isempty(cfg.mcc.method)
    info.method = 'none';
    pAdj = pAlready;
    return;
end

method = lower(string(cfg.mcc.method));
info.method = char(method);

switch method
    case "fwe_permutation"
        fwe     = cfg.mcc.fwe;
        nPerm   = getfield_default(fwe,'nPerm',2000);
        seed    = getfield_default(fwe,'seed',777);
        verbose = getfield_default(fwe,'verbose',true);
        step    = getfield_default(fwe,'progressStep',200);

        [pAdj, info2] = fwe_maxT_unpaired(devData, stdData, nPerm, seed, verbose, step, ttestMode);
        info = info2;
        info.method = sprintf('fwe_permutation(maxT-unpaired, nPerm=%d)', nPerm);

    otherwise
        % Fall back to generic p-map MCC (FDR/Holm/Bonferroni/none)
        [pAdj, info2] = apply_mcc_pmap(pAlready, cfg, tag);
        info = info2;
end
end


function [pFWE, info] = fwe_maxT_unpaired(devData, stdData, nPerm, seed, verbose, progressStep, ttestMode)
% fwe_maxT_unpaired
% Permutation FWE correction using max-|t| for UNPAIRED trial data by label-shuffling.
%
% Inputs:
%   devData/stdData: [chan x time x trials]
%   ttestMode: "student" (pooled var) or "welch"
%
% Output:
%   pFWE: [chan x time] FWE-corrected p-values via max-|t| distribution

rng(seed);

devData = double(devData);
stdData = double(stdData);

n1 = size(devData,3);
n2 = size(stdData,3);

% Observed t-map
[tObs, dfObs] = unpaired_tmap(devData, stdData, ttestMode);
absTobs = abs(tObs(:));

% Build pooled data for permutations
ALL = cat(3, devData, stdData);            % [chan x time x (n1+n2)]
N   = n1 + n2;

maxT = zeros(nPerm,1);

for ip = 1:nPerm
    permIdx = randperm(N);
    grp1 = permIdx(1:n1);
    grp2 = permIdx(n1+1:end);

    X1 = ALL(:,:,grp1);
    X2 = ALL(:,:,grp2);

    tPerm = unpaired_tmap(X1, X2, ttestMode);
    maxT(ip) = max(abs(tPerm(:)));

    if verbose && mod(ip, progressStep) == 0
        fprintf('      FWE maxT (unpaired) perm %d/%d\n', ip, nPerm);
    end
end

% FWE p-values: pFWE(i) = P(max|t_perm| >= |t_obs(i)|)
% Add-one smoothing
maxTsorted = sort(maxT);
pFWE_vec = zeros(size(absTobs));

for i = 1:numel(absTobs)
    nGE = sum(maxTsorted >= absTobs(i));
    pFWE_vec(i) = (1 + nGE) / (1 + nPerm);
end

pFWE = reshape(pFWE_vec, size(tObs));

info = struct();
info.nPerm = nPerm;
info.seed  = seed;
info.df    = dfObs;  % informative; Welch df varies per entry, but we return something sensible
end


function [tMap, dfMap] = unpaired_tmap(devData, stdData, mode)
% unpaired_tmap
% Returns t-stat map for unpaired groups:
%   devData/stdData: [chan x time x trials]
% mode: "student" or "welch"
%
% Output:
%   tMap : [chan x time]
%   dfMap: [chan x time] (Welch varies; Student constant)

devData = double(devData);
stdData = double(stdData);

n1 = size(devData,3);
n2 = size(stdData,3);

m1 = mean(devData, 3);
m2 = mean(stdData, 3);

v1 = var(devData, 0, 3);
v2 = var(stdData, 0, 3);

if strcmpi(mode,"welch")
    se2  = v1./n1 + v2./n2;
    tMap = (m1 - m2) ./ sqrt(se2);

    dfMap = (se2.^2) ./ ( (v1.^2)./(n1^2*(n1-1)) + (v2.^2)./(n2^2*(n2-1)) );
else
    sp2  = ((n1-1).*v1 + (n2-1).*v2) ./ (n1+n2-2);
    se2  = sp2 .* (1/n1 + 1/n2);
    tMap = (m1 - m2) ./ sqrt(se2);

    dfMap = (n1+n2-2) * ones(size(tMap));
end

tMap(~isfinite(tMap)) = 0;
dfMap(~isfinite(dfMap)) = max(n1+n2-2,1);
end


function tf = should_apply_mcc(cfg, level)
% should_apply_mcc
% level: 'individual' or 'group'
% Controlled by cfg.mcc.applyTo = 'none'|'individual'|'group'|'both'

tf = false;
if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'applyTo') || isempty(cfg.mcc.applyTo)
    % default behavior if not specified: apply MCC wherever called
    tf = true;
    return;
end

scope = lower(string(cfg.mcc.applyTo));
level = lower(string(level));

switch scope
    case "none"
        tf = false;
    case "both"
        tf = true;
    case "individual"
        tf = (level == "individual");
    case "group"
        tf = (level == "group");
    otherwise
        warning('Unknown cfg.mcc.applyTo="%s". Using BOTH.', scope);
        tf = true;
end
end


function txt = get_sigbar_label_text(cfg, level)
% Returns a label string for the significance text, depending on MCC.
% Respects cfg.pText: if false -> returns "" (so nothing is printed).

% Gate everything by cfg.pText
if ~(isfield(cfg,'pText') && cfg.pText)
    txt = "";
    return;
end

% Pick alpha to display
alphaTxt = [];
if strcmpi(level,'group') && isfield(cfg,'alpha_grp'), alphaTxt = cfg.alpha_grp; end
if strcmpi(level,'individual') && isfield(cfg,'alpha_ind'), alphaTxt = cfg.alpha_ind; end

% No MCC => plain p
if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'method') || isempty(cfg.mcc.method)
    if ~isempty(alphaTxt), txt = sprintf('p < %.2g', alphaTxt); else, txt = "p"; end
    return;
end

method = lower(string(cfg.mcc.method));

switch method
    case "none"
        if ~isempty(alphaTxt), txt = sprintf('p < %.2g', alphaTxt); else, txt = "p"; end

    case "fdr_bh"
        q = 0.05;
        if isfield(cfg.mcc,'q') && ~isempty(cfg.mcc.q), q = cfg.mcc.q; end
        txt = sprintf('p_{FDR} < %.2g', q);

    case "holm"
        a = 0.05;
        if isfield(cfg.mcc,'alpha') && ~isempty(cfg.mcc.alpha), a = cfg.mcc.alpha; end
        txt = sprintf('p_{Holm} < %.2g', a);

    case "bonferroni"
        a = 0.05;
        if isfield(cfg.mcc,'alpha') && ~isempty(cfg.mcc.alpha), a = cfg.mcc.alpha; end
        txt = sprintf('p_{Bonf} < %.2g', a);

    case "fwe_permutation"
        a = 0.05;
        if isfield(cfg.mcc,'alpha') && ~isempty(cfg.mcc.alpha), a = cfg.mcc.alpha; end
        txt = sprintf('p_{FWE} < %.2g', a);

    otherwise
        if ~isempty(alphaTxt), txt = sprintf('p < %.2g', alphaTxt); else, txt = "p"; end
end
end


function plot_mz4_2x2_with_dual_stats(EEGtmpl, tms, FLS, RLS, FLD, RLD, ...
    pMapLocal, maskTimeLocal, pMapGlobal, maskTimeGlobal, ...
    cfg, thisExp, nSub, outPng, pMapInt, maskTimeInt)

if nargin < 16
    pMapInt = [];
    maskTimeInt = [];
end

chLabels = {EEGtmpl.chanlocs.labels};

barTxt = get_sigbar_label_text(cfg,'group');
if strlength(string(barTxt)) > 0
    barLabelLocal  = "Local effect: "  + barTxt;
    barLabelGlobal = "Global effect: " + barTxt;
else
    barLabelLocal  = "Local effect";
    barLabelGlobal = "Global effect";
end

cBlue   = [0.00 0.45 0.85];
cOrange = [0.90 0.45 0.10];
cGrey   = [0.45 0.45 0.45];

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

tms = double(tms(:)');

for ci = 1:4
    chanName = cfg.erpChans{ci};
    ch = find(strcmpi(chLabels, chanName), 1);
    if isempty(ch), ch = 1; end

    nexttile; hold on; grid on;

    mFLS = double(squeeze(mean(FLS(:,ch,:),1)))';
    mRLS = double(squeeze(mean(RLS(:,ch,:),1)))';
    mFLD = double(squeeze(mean(FLD(:,ch,:),1)))';
    mRLD = double(squeeze(mean(RLD(:,ch,:),1)))';

    semFLS = double(squeeze(std(FLS(:,ch,:),0,1)))' ./ sqrt(size(FLS,1));
    semRLS = double(squeeze(std(RLS(:,ch,:),0,1)))' ./ sqrt(size(RLS,1));
    semFLD = double(squeeze(std(FLD(:,ch,:),0,1)))' ./ sqrt(size(FLD,1));
    semRLD = double(squeeze(std(RLD(:,ch,:),0,1)))' ./ sqrt(size(RLD,1));

    draw_sem_fill(tms, mFLS, semFLS, cBlue,   0.15);
    draw_sem_fill(tms, mRLS, semRLS, cBlue,   0.15);
    draw_sem_fill(tms, mFLD, semFLD, cOrange, 0.15);
    draw_sem_fill(tms, mRLD, semRLD, cOrange, 0.15);

    h1 = plot(tms, mFLS, '--', 'Color', cBlue,   'LineWidth', 1.8);
    h2 = plot(tms, mRLS, '-',  'Color', cBlue,   'LineWidth', 1.8);
    h3 = plot(tms, mFLD, '--', 'Color', cOrange, 'LineWidth', 1.8);
    h4 = plot(tms, mRLD, '-',  'Color', cOrange, 'LineWidth', 1.8);

    xline(0,'k-','LineWidth',1);
    yline(0,'k-','LineWidth',1);

    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    xlim([tms(1) tms(end)]);
    ylim(cfg.ylim);
    title(sprintf('ERP at %s', upper(chanName)), 'FontWeight','bold');

    legend([h1 h2 h3 h4], ...
        {'Frequent local standard','Rare local standard', ...
        'Frequent local deviant','Rare local deviant'}, ...
        'Location','southwest');

    tmsR = tms(:)';

    % ===================== LOCAL EFFECT (TOP BLUE) =====================
    pLocCh = double(pMapLocal(ch,:)); pLocCh = pLocCh(:)';
    mLocR  = logical(maskTimeLocal(:)');

    showMaskLoc = mLocR & (tmsR >= cfg.showFromMs) & (pLocCh <= cfg.alpha_grp);
    showMaskLoc = bridge_small_gaps(showMaskLoc, cfg.barMaxGap);
    showMaskLoc = keep_min_run_length(showMaskLoc, cfg.barMinConsec);

    yBarTop = cfg.ylim(2) - cfg.barYoffset_uV;

    if any(showMaskLoc)
        d = diff([false showMaskLoc false]);
        rs = find(d==1);
        re = find(d==-1)-1;

        for rr = 1:numel(rs)
            idx = rs(rr):re(rr);
            plot(tmsR(idx), yBarTop*ones(size(idx)), '-', ...
                'Color',[0 0 1], 'LineWidth', 6, 'HandleVisibility','off');
        end

        if ~isempty(rs)
            text(tmsR(rs(1)), yBarTop+0.4, barLabelLocal, ...
                'Color',[0 0 1], 'FontSize', 9, 'FontWeight','bold', ...
                'VerticalAlignment','bottom');
        end
    end

    % ===== TRUE OMNIBUS INTERACTION (GREY, ONLY TOP/LOCAL SIDE, BELOW TOP BLUE) =====
    if ~isempty(pMapInt) && ~isempty(maskTimeInt)
        pIntCh = double(pMapInt(ch,:)); pIntCh = pIntCh(:)';
        mIntR  = logical(maskTimeInt(:)');

        showMaskInt = mIntR & (tmsR >= cfg.showFromMs) & (pIntCh <= cfg.alpha_grp);
        showMaskInt = bridge_small_gaps(showMaskInt, cfg.barMaxGap);
        showMaskInt = keep_min_run_length(showMaskInt, cfg.barMinConsec);

        if any(showMaskInt)
            yBarTopInt = yBarTop - 0.55;   % BELOW top blue bar

            d = diff([false showMaskInt false]);
            rs = find(d==1);
            re = find(d==-1)-1;

            for rr = 1:numel(rs)
                idx = rs(rr):re(rr);
                plot(tmsR(idx), yBarTopInt*ones(size(idx)), '-', ...
                    'Color', cGrey, 'LineWidth', 5, 'HandleVisibility','off');
            end
        end
    end

    % ===================== GLOBAL EFFECT (BOTTOM BLUE) =====================
    pGloCh = double(pMapGlobal(ch,:)); pGloCh = pGloCh(:)';
    mGloR  = logical(maskTimeGlobal(:)');

    showMaskGlo = mGloR & (tmsR >= cfg.showFromMs) & (pGloCh <= cfg.alpha_grp);
    showMaskGlo = bridge_small_gaps(showMaskGlo, cfg.barMaxGap);
    showMaskGlo = keep_min_run_length(showMaskGlo, cfg.barMinConsec);

    yBarBot = cfg.ylim(1) + cfg.barYoffset_uV;

    if any(showMaskGlo)
        d = diff([false showMaskGlo false]);
        rs = find(d==1);
        re = find(d==-1)-1;

        for rr = 1:numel(rs)
            idx = rs(rr):re(rr);
            plot(tmsR(idx), yBarBot*ones(size(idx)), '-', ...
                'Color',[0 0 1], 'LineWidth', 6, 'HandleVisibility','off');
        end

        if ~isempty(rs)
            text(tmsR(rs(1)), yBarBot-0.4, barLabelGlobal, ...
                'Color',[0 0 1], 'FontSize', 9, 'FontWeight','bold', ...
                'VerticalAlignment','top');
        end
    end

    box off;
end

if strcmpi(thisExp,'E1')
    sgtitle_experiment = 'Counting';
else
    sgtitle_experiment = 'Visual interference';
end

sgtitle(sprintf('Event Related Potentials | %s', sgtitle_experiment), 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end



function draw_sem_fill(x, m, sem, faceColor, faceAlphaVal)
x   = double(x(:)');
m   = double(m(:)');
sem = double(sem(:)');

if isempty(x) || isempty(m) || isempty(sem), return; end
if ~(numel(x)==numel(m) && numel(m)==numel(sem))
    warning('draw_sem_fill: size mismatch. Skipping fill.');
    return;
end

ylow  = m - sem;
yhigh = m + sem;

fill([x fliplr(x)], [ylow fliplr(yhigh)], faceColor, ...
    'FaceAlpha', faceAlphaVal, 'EdgeColor','none', 'HandleVisibility','off');
end


%% NEW ERP plots puplication style

function eeg_export_pupilstyle_2panel_extras(EEGtmpl, tms, ...
    Lstd, Ldev, Gstd, Gdev, ...
    pLoc_plot, mLoc_plot, pGlo_plot, mGlo_plot, ...
    FLSsub, RLSsub, FLDsub, RLDsub, ...
    cfg, thisExp, nSub, nTrialLocalStd, nTrialLocalDev, nTrialGlobalStd, nTrialGlobalDev, ...
    outDir, pInt_plot, mInt_plot)

% eeg_export_pupilstyle_2panel_extras
% Creates EXTRA (non-invasive) EEG figures in the SAME fixed-size style as your pupil 2-panel plots:
%   - Bekinschtein-style (Std vs Dev) Local + Global, split into 2-channel figures
%   - Mazancieux-style (4 conditions) with Local(top) + Global(bottom) bars, split into 2-channel figures
%
% Does NOT modify your existing outputs/pipeline; just exports additional PNGs.

% ----- fixed canvas (same as your pupil dual) -----
FIG_W_IN = 3735/300;   % 12.45 in
FIG_H_IN = 1247/300;   % 4.1567 in
DPI      = 300;

% Channel pairs from cfg.erpChans (assumes 4 entries: {'Fz','FC1','Cz','Pz'})
chList = cfg.erpChans;
if numel(chList) < 4
    error('cfg.erpChans must have at least 4 channels (e.g., {''Fz'',''FC1'',''Cz'',''Pz''}).');
end
pairs = { chList(1:2), chList(3:4) };

% Friendly experiment string (matches your titles)
if strcmpi(thisExp,'E1')
    expName = 'Counting';
else
    expName = 'Visual interference';
end

% Ensure row vector time
tms = double(tms(:)');

% ===== Bekinschtein-style: Local and Global (2 conditions) =====
for pi = 1:numel(pairs)
    pair = pairs{pi};
    tag  = sprintf('%s_%s', pair{1}, pair{2});

    % Local (Std vs Dev)
    outPng = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_ERP_LOCAL_2panel_%s.png', thisExp, tag));
    eeg_plot_erps_1x2_fixed(EEGtmpl, tms, Lstd, Ldev, pLoc_plot, mLoc_plot, cfg, ...
        expName, 'Local', nSub, nTrialLocalStd, nTrialLocalDev, pair, outPng, FIG_W_IN, FIG_H_IN, DPI, ...
        pInt_plot, mInt_plot);

    % Global (Std vs Dev)
    outPng = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_ERP_GLOBAL_2panel_%s.png', thisExp, tag));
    eeg_plot_erps_1x2_fixed(EEGtmpl, tms, Gstd, Gdev, pGlo_plot, mGlo_plot, cfg, ...
        expName, 'Global', nSub, nTrialGlobalStd, nTrialGlobalDev, pair, outPng, FIG_W_IN, FIG_H_IN, DPI, ...
        pInt_plot, mInt_plot);

end

% ===== Mazancieux-style: 4 conditions + Local(top)/Global(bottom) =====
for pi = 1:numel(pairs)
    pair = pairs{pi};
    tag  = sprintf('%s_%s', pair{1}, pair{2});

    outPng = fullfile(outDir, sprintf('Mazancieux_STATS_%s_ERP4_LOCAL_GLOBAL_2panel_%s.png', thisExp, tag));
    eeg_plot_mz4_1x2_fixed(EEGtmpl, tms, FLSsub, RLSsub, FLDsub, RLDsub, ...
        pLoc_plot, mLoc_plot, pGlo_plot, mGlo_plot, ...
        cfg, expName, nSub, pair, outPng, FIG_W_IN, FIG_H_IN, DPI, ...
        pInt_plot, mInt_plot);
end
end


% =====================================================================
%  Bekinschtein-style 1x2 (two channels) fixed-size exporter
% =====================================================================
function eeg_plot_erps_1x2_fixed(EEGtmpl, tms, STD, DEV, pMap, maskTime, cfg, ...
    expName, whichEffect, nSub, nStdTrials, nDevTrials, chanPair, outPng, FIG_W_IN, FIG_H_IN, DPI, pMapInt, maskTimeInt)

if nargin < 18
    pMapInt = [];
    maskTimeInt = [];
end

chLabels = {EEGtmpl.chanlocs.labels};

fig = figure('Color','w', 'Units','inches', 'Position',[0.5 0.5 FIG_W_IN FIG_H_IN], ...
    'InvertHardcopy','off');
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 FIG_W_IN FIG_H_IN]);
set(fig,'PaperSize',[FIG_W_IN FIG_H_IN]);
set(fig,'PaperPositionMode','manual');

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
tl.Units = 'normalized';
tl.OuterPosition = [0 0 1 1];

for iAx = 1:2
    chanName = chanPair{iAx};
    ch = find(strcmpi(chLabels, chanName), 1);
    if isempty(ch), ch = 1; end

    ax = nexttile; hold(ax,'on'); grid(ax,'on');
    set(ax,'FontSize',11,'Color','w');

    % Mean + SEM across subjects
    mStd  = squeeze(mean(STD(:,ch,:),1));
    mDev  = squeeze(mean(DEV(:,ch,:),1));
    semStd = eeg_sem(STD(:,ch,:));
    semDev = eeg_sem(DEV(:,ch,:));

    % Shaded SEM
    eeg_draw_sem_fill(ax, tms, mStd, semStd, [0 0.45 0.85], 0.15);
    eeg_draw_sem_fill(ax, tms, mDev, semDev, [0.90 0.45 0.10], 0.15);

    % Lines
    plot(ax, tms, mStd, 'LineWidth', 2, 'Color', [0 0.45 0.85]);
    plot(ax, tms, mDev, 'LineWidth', 2, 'Color', [0.90 0.45 0.10]);

    xline(ax,0,'k-','LineWidth',1);
    yline(ax,0,'k-','LineWidth',1);

    xlabel(ax,'Time (ms)');
    ylabel(ax,'Amplitude (\muV)');
    xlim(ax,[tms(1) tms(end)]);
    ylim(ax,cfg.ylim);

    title(ax, sprintf('ERP at %s', upper(chanName)), 'FontWeight','bold');

    % Legend labels depend on effect
    if strcmpi(whichEffect,'Local')
        legStd = 'Local standard'; legDev = 'Local deviant';
    else
        legStd = 'Global standard'; legDev = 'Global deviant';
    end
    legend(ax,{legStd, legDev}, 'Location','southwest');

    % ===== sig bar (time-accurate): maskTime AND pMap(chan,:) =====
    tmsR  = tms(:)';
    maskR = logical(maskTime(:)');

    pch = double(pMap(ch,:)); pch = pch(:)';
    showMask = maskR & (tmsR >= cfg.showFromMs) & (pch <= cfg.alpha_grp);
    showMask = eeg_bridge_small_gaps(showMask, cfg.barMaxGap);
    showMask = eeg_keep_min_run_length(showMask, cfg.barMinConsec);

    if any(showMask)
        yBar = cfg.ylim(2) - cfg.barYoffset_uV;

        [rs,re] = eeg_sig_segments(showMask);
        for rr = 1:numel(rs)
            idx = rs(rr):re(rr);
            plot(ax, tmsR(idx), yBar*ones(size(idx)), '-', 'LineWidth', 6, 'Color', [0 0 1], 'HandleVisibility','off');
        end

        % label BELOW the bar (avoids subplot title)
        if ~isempty(rs)
            if cfg.pText == true
                text(ax, tmsR(rs(1)), yBar-0.35, sprintf('%s effect', whichEffect), ...
                    'Color',[0 0 1],'FontSize',9,'FontWeight','bold', ...
                    'VerticalAlignment','top','HorizontalAlignment','left','Clipping','on');
            end
        end
    end
    % ===== grey interaction bar =====
    if ~isempty(pMapInt) && ~isempty(maskTimeInt)
        pchInt = double(pMapInt(ch,:)); pchInt = pchInt(:)';
        mIntR  = logical(maskTimeInt(:)');

        showMaskInt = mIntR & (tmsR >= cfg.showFromMs) & (pchInt <= cfg.alpha_grp);
        showMaskInt = eeg_bridge_small_gaps(showMaskInt, cfg.barMaxGap);
        showMaskInt = eeg_keep_min_run_length(showMaskInt, cfg.barMinConsec);

        if any(showMaskInt)
            yBarInt = cfg.ylim(2) - cfg.barYoffset_uV - 0.55;

            [rsInt,reInt] = eeg_sig_segments(showMaskInt);
            for rr = 1:numel(rsInt)
                idx = rsInt(rr):reInt(rr);
                plot(ax, tmsR(idx), yBarInt*ones(size(idx)), '-', ...
                    'LineWidth', 5, 'Color', [0.45 0.45 0.45], 'HandleVisibility','off');
            end
        end
    end
end

sgtitle(sprintf('Event Related Potentials | %s | %s effects', expName, whichEffect), 'FontWeight','bold');

% fixed-size export
print(fig, outPng, '-dpng', sprintf('-r%d', DPI));
close(fig); drawnow;
end


% =====================================================================
%  Mazancieux-style 4-condition 1x2 fixed-size exporter (dual bars)
% =====================================================================
function eeg_plot_mz4_1x2_fixed(EEGtmpl, tms, FLS, RLS, FLD, RLD, ...
    pMapLocal, maskTimeLocal, pMapGlobal, maskTimeGlobal, ...
    cfg, expName, nSub, chanPair, outPng, FIG_W_IN, FIG_H_IN, DPI, ...
    pMapInt, maskTimeInt)

if nargin < 20
    pMapInt = [];
    maskTimeInt = [];
end

chLabels = {EEGtmpl.chanlocs.labels};

fig = figure('Color','w', 'Units','inches', 'Position',[0.5 0.5 FIG_W_IN FIG_H_IN], ...
    'InvertHardcopy','off');
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 FIG_W_IN FIG_H_IN]);
set(fig,'PaperSize',[FIG_W_IN FIG_H_IN]);
set(fig,'PaperPositionMode','manual');

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
tl.Units = 'normalized';
tl.OuterPosition = [0 0 1 1];

cBlue   = [0.00 0.45 0.85];
cOrange = [0.90 0.45 0.10];
cGrey   = [0.45 0.45 0.45];

for iAx = 1:2
    chanName = chanPair{iAx};
    ch = find(strcmpi(chLabels, chanName), 1);
    if isempty(ch), ch = 1; end

    ax = nexttile; hold(ax,'on'); grid(ax,'on');
    set(ax,'FontSize',11,'Color','w');

    mFLS = double(squeeze(mean(FLS(:,ch,:),1)))';
    mRLS = double(squeeze(mean(RLS(:,ch,:),1)))';
    mFLD = double(squeeze(mean(FLD(:,ch,:),1)))';
    mRLD = double(squeeze(mean(RLD(:,ch,:),1)))';

    semFLS = eeg_sem(FLS(:,ch,:));
    semRLS = eeg_sem(RLS(:,ch,:));
    semFLD = eeg_sem(FLD(:,ch,:));
    semRLD = eeg_sem(RLD(:,ch,:));

    eeg_draw_sem_fill(ax, tms, mFLS, semFLS, cBlue,   0.15);
    eeg_draw_sem_fill(ax, tms, mRLS, semRLS, cBlue,   0.15);
    eeg_draw_sem_fill(ax, tms, mFLD, semFLD, cOrange, 0.15);
    eeg_draw_sem_fill(ax, tms, mRLD, semRLD, cOrange, 0.15);

    h1 = plot(ax, tms, mFLS, '--', 'Color', cBlue,   'LineWidth', 1.8);
    h2 = plot(ax, tms, mRLS, '-',  'Color', cBlue,   'LineWidth', 1.8);
    h3 = plot(ax, tms, mFLD, '--', 'Color', cOrange, 'LineWidth', 1.8);
    h4 = plot(ax, tms, mRLD, '-',  'Color', cOrange, 'LineWidth', 1.8);

    xline(ax,0,'k-','LineWidth',1);
    yline(ax,0,'k-','LineWidth',1);

    xlabel(ax,'Time (ms)');
    ylabel(ax,'Amplitude (\muV)');
    xlim(ax,[tms(1) tms(end)]);
    ylim(ax,cfg.ylim);
    title(ax, sprintf('ERP at %s', upper(chanName)), 'FontWeight','bold');

    legend(ax,[h1 h2 h3 h4], ...
        {'Frequent local standard','Rare local standard', ...
        'Frequent local deviant','Rare local deviant'}, ...
        'Location','southwest');

    tmsR = tms(:)';

    % ===================== LOCAL EFFECT (TOP BLUE) =====================
    pLocCh = double(pMapLocal(ch,:)); pLocCh = pLocCh(:)';
    mLocR  = logical(maskTimeLocal(:)');

    showLoc = mLocR & (tmsR >= cfg.showFromMs) & (pLocCh <= cfg.alpha_grp);
    showLoc = eeg_bridge_small_gaps(showLoc, cfg.barMaxGap);
    showLoc = eeg_keep_min_run_length(showLoc, cfg.barMinConsec);

    yTop = cfg.ylim(2) - cfg.barYoffset_uV;

    if any(showLoc)
        [rs,re] = eeg_sig_segments(showLoc);
        for rr=1:numel(rs)
            idx = rs(rr):re(rr);
            plot(ax, tmsR(idx), yTop*ones(size(idx)), '-', ...
                'Color',[0 0 1], 'LineWidth', 6, 'HandleVisibility','off');
        end

        if ~isempty(rs)
            text(ax, tmsR(rs(1)), yTop-0.35, 'Local effect', ...
                'Color',[0 0 1],'FontSize',9,'FontWeight','bold', ...
                'VerticalAlignment','top','HorizontalAlignment','left','Clipping','on');
        end
    end

    % ===== TRUE OMNIBUS INTERACTION (GREY, ONLY TOP/LOCAL SIDE, ABOVE TOP BLUE) =====
    if ~isempty(pMapInt) && ~isempty(maskTimeInt)
        pIntCh = double(pMapInt(ch,:)); pIntCh = pIntCh(:)';
        mIntR  = logical(maskTimeInt(:)');

        showLocInt = mIntR & (tmsR >= cfg.showFromMs) & (pIntCh <= cfg.alpha_grp);
        showLocInt = eeg_bridge_small_gaps(showLocInt, cfg.barMaxGap);
        showLocInt = eeg_keep_min_run_length(showLocInt, cfg.barMinConsec);

        if any(showLocInt)
            yTopInt = yTop + 0.55;   % ABOVE top blue bar
            [rs,re] = eeg_sig_segments(showLocInt);
            for rr=1:numel(rs)
                idx = rs(rr):re(rr);
                plot(ax, tmsR(idx), yTopInt*ones(size(idx)), '-', ...
                    'Color', cGrey, 'LineWidth', 5, 'HandleVisibility','off');
            end
        end
    end

    % ===================== GLOBAL EFFECT (BOTTOM BLUE) =====================
    pGloCh = double(pMapGlobal(ch,:)); pGloCh = pGloCh(:)';
    mGloR  = logical(maskTimeGlobal(:)');

    showGlo = mGloR & (tmsR >= cfg.showFromMs) & (pGloCh <= cfg.alpha_grp);
    showGlo = eeg_bridge_small_gaps(showGlo, cfg.barMaxGap);
    showGlo = eeg_keep_min_run_length(showGlo, cfg.barMinConsec);

    yBot = cfg.ylim(1) + cfg.barYoffset_uV;

    if any(showGlo)
        [rs,re] = eeg_sig_segments(showGlo);
        for rr=1:numel(rs)
            idx = rs(rr):re(rr);
            plot(ax, tmsR(idx), yBot*ones(size(idx)), '-', ...
                'Color',[0 0 1], 'LineWidth', 6, 'HandleVisibility','off');
        end

        if ~isempty(rs)
            text(ax, tmsR(rs(1)), yBot+0.35, 'Global effect', ...
                'Color',[0 0 1],'FontSize',9,'FontWeight','bold', ...
                'VerticalAlignment','bottom','HorizontalAlignment','left','Clipping','on');
        end
    end

    box(ax,'off');
end

sgtitle(sprintf('Event Related Potentials | %s', expName), 'FontWeight','bold');

print(fig, outPng, '-dpng', sprintf('-r%d', DPI));
close(fig); drawnow;
end


% =====================================================================
%  Small helpers (self-contained; won't interfere with your code)
% =====================================================================
function sem = eeg_sem(Xsub)
% Xsub expected size: [sub x 1 x time] or [sub x time]
X = double(squeeze(Xsub)); % -> [sub x time]
n = size(X,1);
if n < 2
    sem = zeros(1,size(X,2));
    return;
end
% robust std (omitnan if available)
try
    s = std(X,0,1,'omitnan');
catch
    X(~isfinite(X)) = NaN;
    s = nanstd(X,0,1);
end
sem = s ./ sqrt(n);
sem = sem(:)'; % row
end

function eeg_draw_sem_fill(ax, x, m, sem, faceColor, faceAlphaVal)
x   = double(x(:)');
m   = double(m(:)');
sem = double(sem(:)');

ylow  = m - sem;
yhigh = m + sem;

fill(ax, [x fliplr(x)], [ylow fliplr(yhigh)], faceColor, ...
    'FaceAlpha', faceAlphaVal, 'EdgeColor','none', 'HandleVisibility','off');
end

function mask = eeg_keep_min_run_length(mask, minLen)
mask = logical(mask(:)');
if isempty(mask) || minLen <= 1, return; end
[rs,re] = eeg_sig_segments(mask);
out = false(size(mask));
for i=1:numel(rs)
    if (re(i)-rs(i)+1) >= minLen
        out(rs(i):re(i)) = true;
    end
end
mask = out;
end

function mask = eeg_bridge_small_gaps(mask, maxGap)
mask = logical(mask(:)');
if isempty(mask) || maxGap <= 0, return; end
x = mask;
[rs,re] = eeg_sig_segments(x);
if numel(rs) < 2, mask = x; return; end

for i=1:numel(rs)-1
    gapStart = re(i)+1;
    gapEnd   = rs(i+1)-1;
    gapLen   = gapEnd-gapStart+1;
    if gapLen > 0 && gapLen <= maxGap
        x(gapStart:gapEnd) = true;
    end
end
mask = x;
end

function [starts, ends] = eeg_sig_segments(mask)
mask = logical(mask(:)');
d = diff([false mask false]);
starts = find(d==1);
ends   = find(d==-1)-1;
end

%% Topos per condition

function eeg_export_topo_snaps_conditions(EEGtmpl, tms, Lstd, Ldev, Gstd, Gdev, cfg, thisExp, outDir)
% eeg_export_topo_snaps_conditions
% Exports 4 topography-evolution figures (one per condition):
%   Local Standard, Local Deviant, Global Standard, Global Deviant
% Each figure is a 1xN grid of topoplots at cfg.snapTimes_ms.
%
% Inputs can be either:
%   - [sub x chan x time]  (recommended; we average across subjects)
%   - [chan x time]        (already averaged)

if ~exist(outDir,'dir'), mkdir(outDir); end

% Convert to [chan x time] maps (averaged across subjects if needed)
M_LS = eeg_as_chan_time(Lstd);
M_LD = eeg_as_chan_time(Ldev);
M_GS = eeg_as_chan_time(Gstd);
M_GD = eeg_as_chan_time(Gdev);

% Choose color limits for CONDITIONS:
% Prefer cfg.caxis_erpTopo_cond if you add it; otherwise use cfg.caxis_erpTopo;
% otherwise auto from the 4 condition maps at snap times.
clim = [];
if isfield(cfg,'caxis_erpTopo_cond') && ~isempty(cfg.caxis_erpTopo_cond)
    clim = cfg.caxis_erpTopo_cond;
elseif isfield(cfg,'caxis_erpTopo') && ~isempty(cfg.caxis_erpTopo)
    clim = cfg.caxis_erpTopo;
else
    % Auto symmetric clim from the 4 condition maps (at snap times)
    snap = cfg.snapTimes_ms(:)';
    vmax = 0;
    for s = snap
        [~, ti] = min(abs(tms - s));
        vmax = max(vmax, max(abs([M_LS(:,ti); M_LD(:,ti); M_GS(:,ti); M_GD(:,ti)])));
    end
    if ~isfinite(vmax) || vmax <= 0, vmax = 1; end
    clim = [-vmax vmax];
end

% Output names + titles
items = { ...
    M_LS, sprintf('%s | ERP scalp maps | Local Standard',  thisExp), fullfile(outDir, sprintf('COND_TOPO_%s_LOCAL_STANDARD.png',  thisExp)); ...
    M_LD, sprintf('%s | ERP scalp maps | Local Deviant',   thisExp), fullfile(outDir, sprintf('COND_TOPO_%s_LOCAL_DEVIANT.png',   thisExp)); ...
    M_GS, sprintf('%s | ERP scalp maps | Global Standard', thisExp), fullfile(outDir, sprintf('COND_TOPO_%s_GLOBAL_STANDARD.png', thisExp)); ...
    M_GD, sprintf('%s | ERP scalp maps | Global Deviant',  thisExp), fullfile(outDir, sprintf('COND_TOPO_%s_GLOBAL_DEVIANT.png',  thisExp))  ...
    };

for i = 1:size(items,1)
    mapCT    = items{i,1};
    bigTitle = items{i,2};
    outPng   = items{i,3};

    eeg_plot_topo_snaps_condition(EEGtmpl, tms, mapCT, cfg, bigTitle, outPng, clim);
end
end


function eeg_plot_topo_snaps_condition(EEGtmpl, tms, mapCT, cfg, bigTitle, outPng, clim)
% Compact one-row topoplot export with minimal horizontal gaps

snap   = cfg.snapTimes_ms(:)';
chlocs = EEGtmpl.chanlocs;
nSnap  = numel(snap);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.38]);

% ---------- manual compact layout ----------
leftMargin   = 0.012;
rightMargin  = 0.012;
topMargin    = 0.16;
bottomMargin = 0.20;   % leave room for colorbar
gapX         = 0;  % smaller = tighter
usableW = 1 - leftMargin - rightMargin - (nSnap-1)*gapX;
axW = usableW / nSnap;
axH = 0.58;            % large heads
y0  = bottomMargin + 0.08;

ax = gobjects(nSnap,1);

for i = 1:nSnap
    x0 = leftMargin + (i-1)*(axW + gapX);

    ax(i) = axes('Parent',fig,'Units','normalized', ...
        'Position',[x0 y0 axW axH]);

    [~, ti] = min(abs(tms - snap(i)));

    topoplot(double(mapCT(:,ti)), chlocs, ...
        'electrodes','on', ...
        'numcontour', cfg.topoContours, ...
        'gridscale',  cfg.topoGrid, ...
        'shading',    cfg.topoShading, ...
        'whitebk',    'off');

    if ~isempty(clim)
        caxis(ax(i), clim);
    end

    title(ax(i), sprintf('%d ms', snap(i)), ...
        'FontSize', 12, 'FontWeight','bold');
end

sgtitle(fig, bigTitle, 'FontWeight','bold');

cb = colorbar('southoutside');
cb.Label.String = '\muV';
cb.Label.FontSize = 12;
cb.Position = [0.37 0.11 0.25 0.045];

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end


function M = eeg_as_chan_time(X)
% Accepts [sub x chan x time] OR [chan x time] and returns [chan x time]
X = double(X);
if ndims(X) == 3
    M = squeeze(mean(X, 1));  % mean across subjects -> [chan x time]
elseif ndims(X) == 2
    M = X;
else
    error('eeg_as_chan_time: expected [sub x chan x time] or [chan x time].');
end
end


function eeg_expand_tiled_axes(ax, left, bottom, width, height)
% Same behavior as your expand_tiled_axes, but uniquely named.
ax = ax(isgraphics(ax));
if isempty(ax), return; end

pos = cell2mat(get(ax,'Position'));
minX = min(pos(:,1));
minY = min(pos(:,2));
maxX = max(pos(:,1)+pos(:,3));
maxY = max(pos(:,2)+pos(:,4));

curW = maxX - minX;
curH = maxY - minY;

for i = 1:size(pos,1)
    x = (pos(i,1)-minX)/curW;
    y = (pos(i,2)-minY)/curH;
    w = pos(i,3)/curW;
    h = pos(i,4)/curH;

    newPos = [left + x*width, bottom + y*height, w*width, h*height];
    set(ax(i),'Position',newPos);
end
end

%% ========================================================================
%  NEW 2-panel summary exports (EEG) in the SAME style as your PUPIL exports
% ========================================================================
% EEG summary figures in EXACT pupil-style canvas/layout
% ------------------------------------------------------------------------
% Creates 9 fixed-size summary figures (3735x1247 px @300 dpi):
%   1) Bekinschtein_STATS_SSI_MATRIX_2panel_E1.png
%   2) Bekinschtein_STATS_SSI_MATRIX_2panel_E2.png
%   3) Bekinschtein_STATS_GROUP_PVAL_MATRIX_2panel_E1.png
%   4) Bekinschtein_STATS_GROUP_PVAL_MATRIX_2panel_E2.png
%   5) Bekinschtein_STATS_BinnedPvals_SubjectMatrix_2panel_E1.png
%   6) Bekinschtein_STATS_BinnedPvals_SubjectMatrix_2panel_E2.png
%   7) Bekinschtein_STATS_BinnedPvals_ElectrodeMatrix_2panel_E1.png
%   8) Bekinschtein_STATS_BinnedPvals_ElectrodeMatrix_2panel_E2.png
%   9) Bekinschtein_STATS_SSI_MEAN_INDEX_2panel_E1_E2.png
%
% Matches the pupil export style:
%   - same canvas
%   - same tiledlayout spacing
%   - same colorbar placement
%   - same text positioning
%   - same font sizing logic
%   - same "zoomed-out" feel
% ========================================================================

function eeg_export_9_extra_zoomed_figures(outroot, grandFold, SSI, GROUPPMAP, cfg, cfgBIN_EEG, excludeSubs, targets)

dpi = 300;
FIG_W_IN = 3735/dpi;
FIG_H_IN = 1247/dpi;

if nargin < 7 || isempty(excludeSubs), excludeSubs = {}; end
if nargin < 8 || isempty(targets)
    targets = ["Frequent Local Standard","Rare Local Standard", ...
        "Frequent Local Deviant","Rare Local Deviant"];
end

if ~exist(grandFold,'dir'), error('grandFold not found: %s', grandFold); end
if ~exist(outroot,'dir'),   error('outroot not found: %s', outroot); end

%% ================================================================
% 1-2) SUBJECT-WISE SSI DARK MATRICES
% ================================================================
pb    = cfg.pBins_matrix(:)';
nBins = numel(pb);
cmapDark = eegzoom_pbin_cmap(cfg, nBins+1);

hideAbove = inf;
if isfield(cfg,'displayHideBinsAboveP') && isfinite(cfg.displayHideBinsAboveP)
    hideAbove = cfg.displayHideBinsAboveP;
end

% ---------- E1 subject SSI ----------
if isfield(SSI,'E1') && ~isempty(SSI.E1.local)
    tms = SSI.E1.times(:)';
    keepPost0 = eegzoom_keep_post0(tms, cfg);
    tms2 = tms(keepPost0);

    E1L = SSI.E1.local(:,keepPost0);
    E1G = SSI.E1.global(:,keepPost0);
    [E1L,E1G] = eegzoom_apply_ssi_hide(E1L,E1G,pb,hideAbove);

    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_imagesc_matrix(ax1, tms2, E1L, ...
        'Event Related Potentials | Counting', 'Local effects', cmapDark, nBins, 'Participants');

    ax2 = nexttile;
    eegzoom_imagesc_matrix(ax2, tms2, E1G, ...
        'Event Related Potentials | Counting', 'Global effects', cmapDark, nBins, 'Participants');

    eegzoom_add_ssi_colorbar(ax2, pb, cmapDark, cfg, [-0.5 nBins+0.5], [0.90 0.21 0.015 0.58]);

    outPng = fullfile(grandFold,'Bekinschtein_STATS_SSI_MATRIX_2panel_E1.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

% ---------- E2 subject SSI ----------
if isfield(SSI,'E2') && ~isempty(SSI.E2.local)
    tms = SSI.E2.times(:)';
    keepPost0 = eegzoom_keep_post0(tms, cfg);
    tms2 = tms(keepPost0);

    E2L = SSI.E2.local(:,keepPost0);
    E2G = SSI.E2.global(:,keepPost0);
    [E2L,E2G] = eegzoom_apply_ssi_hide(E2L,E2G,pb,hideAbove);

    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_imagesc_matrix(ax1, tms2, E2L, ...
        'Event Related Potentials | Visual interference', 'Local effects', cmapDark, nBins, 'Participants');

    ax2 = nexttile;
    eegzoom_imagesc_matrix(ax2, tms2, E2G, ...
        'Event Related Potentials | Visual interference', 'Global effects', cmapDark, nBins, 'Participants');

    eegzoom_add_ssi_colorbar(ax2, pb, cmapDark, cfg, [-0.5 nBins+0.5], [0.90 0.21 0.015 0.58]);

    outPng = fullfile(grandFold,'Bekinschtein_STATS_SSI_MATRIX_2panel_E2.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

%% ================================================================
% 3-4) GROUP DARK MATRICES (electrodes x time)
% ================================================================
if isfield(GROUPPMAP,'E1') && ~isempty(GROUPPMAP.E1.pLoc)
    tms = GROUPPMAP.E1.times(:)';
    keepPost0 = eegzoom_keep_post0(tms, cfg);
    tms2 = tms(keepPost0);

    MLoc = eegzoom_group_dark_matrix(GROUPPMAP.E1.pLoc, GROUPPMAP.E1.mLoc, pb, keepPost0);
    MGlo = eegzoom_group_dark_matrix(GROUPPMAP.E1.pGlo, GROUPPMAP.E1.mGlo, pb, keepPost0);
    [MLoc,MGlo] = eegzoom_apply_ssi_hide(MLoc,MGlo,pb,hideAbove);

    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_imagesc_matrix(ax1, tms2, MLoc, ...
        'Event Related Potentials | Counting', 'Local effects', cmapDark, nBins, 'Electrodes');

    ax2 = nexttile;
    eegzoom_imagesc_matrix(ax2, tms2, MGlo, ...
        'Event Related Potentials | Counting', 'Global effects', cmapDark, nBins, 'Electrodes');

    eegzoom_add_ssi_colorbar(ax2, pb, cmapDark, cfg, [-0.5 nBins+0.5], [0.90 0.21 0.015 0.58]);

    outPng = fullfile(grandFold,'Bekinschtein_STATS_GROUP_PVAL_MATRIX_2panel_E1.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

if isfield(GROUPPMAP,'E2') && ~isempty(GROUPPMAP.E2.pLoc)
    tms = GROUPPMAP.E2.times(:)';
    keepPost0 = eegzoom_keep_post0(tms, cfg);
    tms2 = tms(keepPost0);

    MLoc = eegzoom_group_dark_matrix(GROUPPMAP.E2.pLoc, GROUPPMAP.E2.mLoc, pb, keepPost0);
    MGlo = eegzoom_group_dark_matrix(GROUPPMAP.E2.pGlo, GROUPPMAP.E2.mGlo, pb, keepPost0);
    [MLoc,MGlo] = eegzoom_apply_ssi_hide(MLoc,MGlo,pb,hideAbove);

    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_imagesc_matrix(ax1, tms2, MLoc, ...
        'Event Related Potentials | Visual interference', 'Local effects', cmapDark, nBins, 'Electrodes');

    ax2 = nexttile;
    eegzoom_imagesc_matrix(ax2, tms2, MGlo, ...
        'Event Related Potentials | Visual interference', 'Global effects', cmapDark, nBins, 'Electrodes');

    eegzoom_add_ssi_colorbar(ax2, pb, cmapDark, cfg, [-0.5 nBins+0.5], [0.90 0.21 0.015 0.58]);

    outPng = fullfile(grandFold,'Bekinschtein_STATS_GROUP_PVAL_MATRIX_2panel_E2.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

%% ================================================================
% 5-8) BINNED BLUE-WHITE MATRICES
% ================================================================
BIN = eegzoom_build_binned_matrices(outroot, cfg, cfgBIN_EEG, excludeSubs, targets);
cmapBW = eegzoom_cmap_bluewhite(256);

% ---------- Subject matrix E1 ----------
if isfield(BIN,'E1') && isfield(BIN.E1,'pLocSub') && ~isempty(BIN.E1.pLocSub)
    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_plot_binned_pvals(ax1, BIN.tCenters, BIN.E1.pLocSub, BIN.E1.sigLocSub, ...
        'Event Related Potentials | Counting', 'Local effects', 'Participants', cmapBW, BIN.cMax);

    ax2 = nexttile;
    eegzoom_plot_binned_pvals(ax2, BIN.tCenters, BIN.E1.pGloSub, BIN.E1.sigGloSub, ...
        'Event Related Potentials | Counting', 'Global effects', 'Participants', cmapBW, BIN.cMax);

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

    outPng = fullfile(grandFold,'Bekinschtein_STATS_BinnedPvals_SubjectMatrix_2panel_E1.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

% ---------- Subject matrix E2 ----------
if isfield(BIN,'E2') && isfield(BIN.E2,'pLocSub') && ~isempty(BIN.E2.pLocSub)
    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_plot_binned_pvals(ax1, BIN.tCenters, BIN.E2.pLocSub, BIN.E2.sigLocSub, ...
        'Event Related Potentials | Visual interference', 'Local effects', 'Participants', cmapBW, BIN.cMax);

    ax2 = nexttile;
    eegzoom_plot_binned_pvals(ax2, BIN.tCenters, BIN.E2.pGloSub, BIN.E2.sigGloSub, ...
        'Event Related Potentials | Visual interference', 'Global effects', 'Participants', cmapBW, BIN.cMax);

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

    outPng = fullfile(grandFold,'Bekinschtein_STATS_BinnedPvals_SubjectMatrix_2panel_E2.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

% ---------- Electrode matrix E1 ----------
if isfield(BIN,'E1') && isfield(BIN.E1,'pLocElec') && ~isempty(BIN.E1.pLocElec)
    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_plot_binned_pvals(ax1, BIN.tCenters, BIN.E1.pLocElec, BIN.E1.sigLocElec, ...
        'Event Related Potentials | Counting', 'Local effects', 'Electrodes', cmapBW, BIN.cMax);

    ax2 = nexttile;
    eegzoom_plot_binned_pvals(ax2, BIN.tCenters, BIN.E1.pGloElec, BIN.E1.sigGloElec, ...
        'Event Related Potentials | Counting', 'Global effects', 'Electrodes', cmapBW, BIN.cMax);

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

    outPng = fullfile(grandFold,'Bekinschtein_STATS_BinnedPvals_ElectrodeMatrix_2panel_E1.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

% ---------- Electrode matrix E2 ----------
if isfield(BIN,'E2') && isfield(BIN.E2,'pLocElec') && ~isempty(BIN.E2.pLocElec)
    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    tl.OuterPosition = [0.02 0.08 0.92 0.86];

    ax1 = nexttile;
    eegzoom_plot_binned_pvals(ax1, BIN.tCenters, BIN.E2.pLocElec, BIN.E2.sigLocElec, ...
        'Event Related Potentials | Visual interference', 'Local effects', 'Electrodes', cmapBW, BIN.cMax);

    ax2 = nexttile;
    eegzoom_plot_binned_pvals(ax2, BIN.tCenters, BIN.E2.pGloElec, BIN.E2.sigGloElec, ...
        'Event Related Potentials | Visual interference', 'Global effects', 'Electrodes', cmapBW, BIN.cMax);

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

    outPng = fullfile(grandFold,'Bekinschtein_STATS_BinnedPvals_ElectrodeMatrix_2panel_E2.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

%% ================================================================
% 9) MEAN SSI INDEX LINE PLOT
% ================================================================
if isfield(SSI,'E1') && isfield(SSI,'E2') && ~isempty(SSI.E1.local) && ~isempty(SSI.E2.local)

    tms = SSI.E1.times(:)';
    keepPost0 = eegzoom_keep_post0(tms, cfg);
    tms2 = tms(keepPost0);

    mE1L = mean(SSI.E1.local(:,keepPost0),1,'omitnan');
    mE1G = mean(SSI.E1.global(:,keepPost0),1,'omitnan');
    mE2L = mean(SSI.E2.local(:,keepPost0),1,'omitnan');
    mE2G = mean(SSI.E2.global(:,keepPost0),1,'omitnan');

    yMax = max([mE1L(:); mE1G(:); mE2L(:); mE2G(:)], [], 'omitnan');
    yMax = max(1, 1.05*yMax);

    fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
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

    sgtitle(fig,'Mean individual statistical significance index | Event Related Potentials','FontWeight','bold');

    outPng = fullfile(grandFold,'Bekinschtein_STATS_SSI_MEAN_INDEX_2panel_E1_E2.png');
    eegzoom_print_fixed(fig, outPng, dpi);
    close(fig);
end

fprintf('[eeg_export_9_extra_zoomed_figures] DONE -> %s\n', grandFold);
end

%% ========================================================================
% HELPERS
% ========================================================================

function fig = eegzoom_fig_fixed(FIG_W_IN, FIG_H_IN)
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 FIG_W_IN FIG_H_IN], ...
    'InvertHardcopy','off');
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0 0 FIG_W_IN FIG_H_IN]);
set(fig,'PaperSize',[FIG_W_IN FIG_H_IN]);
set(fig,'PaperPositionMode','manual');
end

function eegzoom_print_fixed(fig, outPng, dpi)
print(fig, outPng, '-dpng', sprintf('-r%d', dpi));
end

function keepPost0 = eegzoom_keep_post0(tms, cfg)
if isfield(cfg,'showFromMs') && ~isempty(cfg.showFromMs)
    keepPost0 = (tms >= cfg.showFromMs);
else
    keepPost0 = (tms >= 0);
end
end

function cmap = eegzoom_pbin_cmap(cfg, nColors)
if isfield(cfg,'binColors') && ~isempty(cfg.binColors) ...
        && size(cfg.binColors,2)==3 && size(cfg.binColors,1)==nColors
    cmap = cfg.binColors;
else
    cmap = jet(max(nColors,2)); %#ok<JET>
end
end

function [A,B] = eegzoom_apply_ssi_hide(A,B,pBins,hideAboveP)
if ~isfinite(hideAboveP), return; end
hideLvl = find(pBins <= hideAboveP, 1, 'last');
if isempty(hideLvl), hideLvl = 0; end
A(A < hideLvl) = 0;
B(B < hideLvl) = 0;
end

function eegzoom_imagesc_matrix(ax, tms, M, titleTop, titleSmall, cmap, nBins, ylabStr)
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

function eegzoom_add_ssi_colorbar(ax, pBins, cmap, cfg, clim, pos)
nBins = numel(pBins);
cb = colorbar(ax);
cb.Position = pos;
colormap(ax, cmap);
cb.Limits = clim;

ticks = 0:nBins;
cb.Ticks = ticks;

tickLabels = strings(size(ticks));
noneLabel = ">0.05";
if isfield(cfg,'noneLabel') && strlength(string(cfg.noneLabel))>0
    noneLabel = string(cfg.noneLabel);
end
tickLabels(1) = noneLabel;
for k = 1:nBins
    tickLabels(k+1) = sprintf('<=%.4g', pBins(k));
end
cb.TickLabels = tickLabels;

lbl = 'p-values';
if isfield(cfg,'binLabel') && strlength(string(cfg.binLabel))>0
    lbl = char(cfg.binLabel);
end
cb.Label.String = lbl;

cb.Label.Rotation = 90;
cb.Label.Units = 'normalized';
cb.Label.Position = [5.8 0.5 0]; % the firt number: smaller-more to left, bigger more to right
cb.Label.HorizontalAlignment = 'center';
cb.Label.VerticalAlignment   = 'middle';
cb.FontSize = 12;
end

function M = eegzoom_group_dark_matrix(pCT, maskTime, pBins, keepPost0)
pCT = pCT(:,keepPost0);
maskTime = maskTime(keepPost0);

M = zeros(size(pCT));
for tt = 1:numel(maskTime)
    if ~maskTime(tt), continue; end
    pv = pCT(:,tt);
    lvl = zeros(size(pv));
    for k = 1:numel(pBins)
        lvl(pv <= pBins(k)) = k;
    end
    M(:,tt) = lvl;
end
end

function eegzoom_plot_binned_pvals(ax, tCenters, P, sigMask, titleTop, titleSmall, ylabStr, cmapBW, cMax)
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

function cmap = eegzoom_cmap_bluewhite(n)
if nargin < 1 || isempty(n), n = 256; end
a = [247,251,255]/255; b = [222,235,247]/255; c = [198,219,239]/255;
d = [158,202,225]/255; e = [107,174,214]/255; f = [66,146,198]/255; g = [33,113,181]/255;
customColors = [g; f; e; d; c; b; a];
cmap = interp1(linspace(0,1,size(customColors,1)), customColors, linspace(0,1,n), 'pchip');
end

%% ========================================================================
% BUILD EEG BINNED MATRICES
% ========================================================================

function BIN = eegzoom_build_binned_matrices(outroot, cfg, cfgBIN, excludeSubs, targets)

BIN = struct();

binSize = cfgBIN.binSizeMS;
binStep = cfgBIN.binStepMS;
winMS   = cfgBIN.winMS;

tStarts  = winMS(1):binStep:(winMS(2)-binSize);
tCenters = tStarts + binSize/2;

BIN.tCenters = tCenters(:)';
if isfield(cfgBIN,'climVals') && numel(cfgBIN.climVals)>=2
    BIN.cMax = cfgBIN.climVals(2);
else
    BIN.cMax = 0.05;
end

if ~isfield(cfgBIN,'alpha_ind') || isempty(cfgBIN.alpha_ind), cfgBIN.alpha_ind = 0.01; end
if ~isfield(cfgBIN,'minElec_ind') || isempty(cfgBIN.minElec_ind), cfgBIN.minElec_ind = 5; end
if ~isfield(cfgBIN,'indTtest') || isempty(cfgBIN.indTtest), cfgBIN.indTtest = "student"; end
if ~isfield(cfgBIN,'doBalance') || isempty(cfgBIN.doBalance), cfgBIN.doBalance = true; end
if ~isfield(cfgBIN,'balanceSeed') || isempty(cfgBIN.balanceSeed), cfgBIN.balanceSeed = 7; end
if ~isfield(cfgBIN,'doFWE') || isempty(cfgBIN.doFWE), cfgBIN.doFWE = true; end
if ~isfield(cfgBIN,'fweMethod') || isempty(cfgBIN.fweMethod), cfgBIN.fweMethod = "holm"; end
if ~isfield(cfgBIN,'alphaFWE') || isempty(cfgBIN.alphaFWE), cfgBIN.alphaFWE = 0.05; end
if ~isfield(cfgBIN,'minConsecBins_ind') || isempty(cfgBIN.minConsecBins_ind), cfgBIN.minConsecBins_ind = 1; end
if ~isfield(cfgBIN,'minConsecBins_grp') || isempty(cfgBIN.minConsecBins_grp), cfgBIN.minConsecBins_grp = 1; end

alphaGrp = 0.05;
if isfield(cfg,'alpha_grp') && ~isempty(cfg.alpha_grp)
    alphaGrp = cfg.alpha_grp;
end

expList = {'E1','E2'};

for eidx = 1:numel(expList)
    expKey = expList{eidx};
    subs = eegzoom_list_subjects(outroot, excludeSubs, expKey);

    if isempty(subs), continue; end

    templateLabels = [];
    templateTimes  = [];

    for i = 1:numel(subs)
        fAll = fullfile(outroot, char(subs(i)), sprintf('auditory_all_LG_%s_%s.set', char(subs(i)), expKey));
        if exist(fAll,'file')
            EEGt = pop_loadset(fAll);
            EEGt = eeg_checkset(EEGt);
            templateLabels = {EEGt.chanlocs.labels};
            templateTimes  = double(EEGt.times(:))';
            break;
        end
    end

    if isempty(templateLabels), continue; end

    nBins = numel(tCenters);
    nChan = numel(templateLabels);

    pLocSub    = nan(numel(subs), nBins);
    pGloSub    = nan(numel(subs), nBins);
    candLocSub = false(numel(subs), nBins);
    candGloSub = false(numel(subs), nBins);

    keepRow = false(numel(subs),1);

    LstdAll = [];
    LdevAll = [];
    GstdAll = [];
    GdevAll = [];

    for i = 1:numel(subs)
        participant = char(subs(i));
        fAll = fullfile(outroot, participant, sprintf('auditory_all_LG_%s_%s.set', participant, expKey));
        if ~exist(fAll,'file'), continue; end

        EEGs = pop_loadset(fAll);
        EEGs = eeg_checkset(EEGs);
        EEGs = eegzoom_reorder_to_template(EEGs, templateLabels);

        times = double(EEGs.times(:))';
        if numel(times) ~= numel(templateTimes) || any(abs(times - templateTimes) > 1e-6)
            warning('Time mismatch in %s %s. Skipping.', participant, expKey);
            continue;
        end

        if ~(isfield(EEGs,'epoch') && isfield(EEGs.epoch,'eventlg_type'))
            continue;
        end

        M = false(4, EEGs.trials);
        for ep = 1:EEGs.trials
            lablist = eegzoom_to_strlist(EEGs.epoch(ep).eventlg_type);
            for k = 1:4
                M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
            end
        end

        idxFLS = find(M(1,:));
        idxRLS = find(M(2,:));
        idxFLD = find(M(3,:));
        idxRLD = find(M(4,:));

        if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
            continue;
        end

        keepRow(i) = true;

        ERP_FLS = mean(double(EEGs.data(:,:,idxFLS)), 3);
        ERP_RLS = mean(double(EEGs.data(:,:,idxRLS)), 3);
        ERP_FLD = mean(double(EEGs.data(:,:,idxFLD)), 3);
        ERP_RLD = mean(double(EEGs.data(:,:,idxRLD)), 3);

        lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3);
        ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3);
        gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3);
        gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3);

        LstdAll(end+1,:,:) = lstd; %#ok<AGROW>
        LdevAll(end+1,:,:) = ldev; %#ok<AGROW>
        GstdAll(end+1,:,:) = gstd; %#ok<AGROW>
        GdevAll(end+1,:,:) = gdev; %#ok<AGROW>

        if cfgBIN.doBalance
            rng(cfgBIN.balanceSeed + 1000*eidx + i);
            [devLoc, stdLoc] = eegzoom_build_balanced_trials(EEGs.data, idxFLD, idxRLD, idxFLS, idxRLS);
            [devGlo, stdGlo] = eegzoom_build_balanced_trials(EEGs.data, idxRLS, idxRLD, idxFLS, idxFLD);
        else
            devLoc = double(EEGs.data(:,:, [idxFLD(:); idxRLD(:)] ));
            stdLoc = double(EEGs.data(:,:, [idxFLS(:); idxRLS(:)] ));
            devGlo = double(EEGs.data(:,:, [idxRLS(:); idxRLD(:)] ));
            stdGlo = double(EEGs.data(:,:, [idxFLS(:); idxFLD(:)] ));
        end

        for b = 1:nBins
            wMask = (times >= tStarts(b)) & (times < (tStarts(b) + binSize));
            if ~any(wMask), continue; end

            XL1 = squeeze(mean(devLoc(:,wMask,:), 2));
            XL2 = squeeze(mean(stdLoc(:,wMask,:), 2));
            XG1 = squeeze(mean(devGlo(:,wMask,:), 2));
            XG2 = squeeze(mean(stdGlo(:,wMask,:), 2));

            pLocCh = eegzoom_unpaired_bin_pvals(XL1, XL2, cfgBIN.indTtest);
            pGloCh = eegzoom_unpaired_bin_pvals(XG1, XG2, cfgBIN.indTtest);

            pLocSub(i,b) = min(pLocCh, [], 'omitnan');
            pGloSub(i,b) = min(pGloCh, [], 'omitnan');

            candLocSub(i,b) = sum(pLocCh <= cfgBIN.alpha_ind) >= cfgBIN.minElec_ind;
            candGloSub(i,b) = sum(pGloCh <= cfgBIN.alpha_ind) >= cfgBIN.minElec_ind;
        end
    end

    subs       = subs(keepRow);
    pLocSub    = pLocSub(keepRow,:);
    pGloSub    = pGloSub(keepRow,:);
    candLocSub = candLocSub(keepRow,:);
    candGloSub = candGloSub(keepRow,:);

    LstdAll = LstdAll(1:sum(keepRow),:,:);
    LdevAll = LdevAll(1:sum(keepRow),:,:);
    GstdAll = GstdAll(1:sum(keepRow),:,:);
    GdevAll = GdevAll(1:sum(keepRow),:,:);

    sigLocSub = false(size(pLocSub));
    sigGloSub = false(size(pGloSub));

    for r = 1:size(pLocSub,1)
        if cfgBIN.doFWE
            sigLocSub(r,:) = eegzoom_fwe_mask(pLocSub(r,:), cfgBIN.alphaFWE, cfgBIN.fweMethod);
            sigGloSub(r,:) = eegzoom_fwe_mask(pGloSub(r,:), cfgBIN.alphaFWE, cfgBIN.fweMethod);
        else
            sigLocSub(r,:) = (pLocSub(r,:) <= cfgBIN.alphaFWE);
            sigGloSub(r,:) = (pGloSub(r,:) <= cfgBIN.alphaFWE);
        end

        sigLocSub(r,:) = sigLocSub(r,:) & candLocSub(r,:);
        sigGloSub(r,:) = sigGloSub(r,:) & candGloSub(r,:);

        sigLocSub(r,:) = eegzoom_keep_min_run_length(sigLocSub(r,:), cfgBIN.minConsecBins_ind);
        sigGloSub(r,:) = eegzoom_keep_min_run_length(sigGloSub(r,:), cfgBIN.minConsecBins_ind);
    end

    pLocElec = nan(nChan, nBins);
    pGloElec = nan(nChan, nBins);

    if ~isempty(LstdAll)
        for b = 1:nBins
            wMask = (templateTimes >= tStarts(b)) & (templateTimes < (tStarts(b) + binSize));
            if ~any(wMask), continue; end

            Astd = squeeze(mean(LstdAll(:,:,wMask), 3));
            Adev = squeeze(mean(LdevAll(:,:,wMask), 3));
            Gstd = squeeze(mean(GstdAll(:,:,wMask), 3));
            Gdev = squeeze(mean(GdevAll(:,:,wMask), 3));

            pLocElec(:,b) = eegzoom_paired_bin_pvals(Adev, Astd);
            pGloElec(:,b) = eegzoom_paired_bin_pvals(Gdev, Gstd);
        end
    end

    sigLocElec = false(size(pLocElec));
    sigGloElec = false(size(pGloElec));

    for ch = 1:size(pLocElec,1)
        if cfgBIN.doFWE
            sigLocElec(ch,:) = eegzoom_fwe_mask(pLocElec(ch,:), cfgBIN.alphaFWE, cfgBIN.fweMethod);
            sigGloElec(ch,:) = eegzoom_fwe_mask(pGloElec(ch,:), cfgBIN.alphaFWE, cfgBIN.fweMethod);
        else
            sigLocElec(ch,:) = (pLocElec(ch,:) <= alphaGrp);
            sigGloElec(ch,:) = (pGloElec(ch,:) <= alphaGrp);
        end

        sigLocElec(ch,:) = eegzoom_keep_min_run_length(sigLocElec(ch,:), cfgBIN.minConsecBins_grp);
        sigGloElec(ch,:) = eegzoom_keep_min_run_length(sigGloElec(ch,:), cfgBIN.minConsecBins_grp);
    end

    BIN.(expKey).subjects   = subs;
    BIN.(expKey).pLocSub    = pLocSub;
    BIN.(expKey).pGloSub    = pGloSub;
    BIN.(expKey).sigLocSub  = sigLocSub;
    BIN.(expKey).sigGloSub  = sigGloSub;

    BIN.(expKey).pLocElec   = pLocElec;
    BIN.(expKey).pGloElec   = pGloElec;
    BIN.(expKey).sigLocElec = sigLocElec;
    BIN.(expKey).sigGloElec = sigGloElec;
    BIN.(expKey).chanLabels = string(templateLabels);
end
end

function subs = eegzoom_list_subjects(outroot, excludeSubs, expKey)
d = dir(fullfile(outroot,'S*'));
d = d([d.isdir]);
names = string({d.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+', 'once'));
subs = names(keep);

if ~isempty(excludeSubs)
    subs = subs(~ismember(lower(subs), lower(string(excludeSubs))));
end

hasFile = false(size(subs));
for i = 1:numel(subs)
    fAll = fullfile(outroot, char(subs(i)), sprintf('auditory_all_LG_%s_%s.set', char(subs(i)), expKey));
    hasFile(i) = exist(fAll,'file') == 2;
end
subs = subs(hasFile);

nums = double(extractAfter(subs,"S"));
[~,ix] = sort(nums);
subs = subs(ix);
end

function EEGout = eegzoom_reorder_to_template(EEGin, templateLabels)
curLabels = {EEGin.chanlocs.labels};
[~, idx] = ismember(templateLabels, curLabels);
if any(idx==0)
    missing = templateLabels(idx==0);
    error('Missing channels: %s', strjoin(missing, ', '));
end
EEGout = EEGin;
EEGout.data     = EEGin.data(idx,:,:);
EEGout.chanlocs = EEGin.chanlocs(idx);
EEGout.nbchan   = numel(idx);
EEGout = eeg_checkset(EEGout);
end

function S = eegzoom_to_strlist(v)
if isstring(v)
    S = v;
elseif ischar(v)
    S = string(v);
elseif iscell(v)
    S = strings(0,1);
    for i=1:numel(v)
        if isstring(v{i}) || ischar(v{i})
            S(end+1,1) = string(v{i}); %#ok<AGROW>
        end
    end
else
    S = strings(0,1);
end
S = strtrim(regexprep(S,'\s+',' '));
end

function [devData, stdData] = eegzoom_build_balanced_trials(DATA, devA, devB, stdA, stdB)
devA = devA(:)'; devB = devB(:)'; stdA = stdA(:)'; stdB = stdB(:)';

nDevMin = min(numel(devA), numel(devB));
nStdMin = min(numel(stdA), numel(stdB));

devA = devA(randperm(numel(devA), nDevMin));
devB = devB(randperm(numel(devB), nDevMin));
stdA = stdA(randperm(numel(stdA), nStdMin));
stdB = stdB(randperm(numel(stdB), nStdMin));

devIdx = [devA devB];
stdIdx = [stdA stdB];

nMin = min(numel(devIdx), numel(stdIdx));
devIdx = devIdx(randperm(numel(devIdx), nMin));
stdIdx = stdIdx(randperm(numel(stdIdx), nMin));

devData = double(DATA(:,:,devIdx));
stdData = double(DATA(:,:,stdIdx));
end

function p = eegzoom_unpaired_bin_pvals(X1, X2, mode)
X1 = double(X1);
X2 = double(X2);

if isvector(X1), X1 = X1(:); end
if isvector(X2), X2 = X2(:); end
if size(X1,1) == 1 && size(X2,1) ~= 1, X1 = X1'; end
if size(X2,1) == 1 && size(X1,1) ~= 1, X2 = X2'; end

if size(X1,2) < 2 || size(X2,2) < 2
    p = ones(size(X1,1),1);
    return;
end

n1 = size(X1,2);
n2 = size(X2,2);

m1 = mean(X1, 2, 'omitnan');
m2 = mean(X2, 2, 'omitnan');

v1 = var(X1, 0, 2, 'omitnan');
v2 = var(X2, 0, 2, 'omitnan');

if strcmpi(mode,"welch")
    se2 = v1./n1 + v2./n2;
    t   = (m1 - m2) ./ sqrt(se2);
    df  = (se2.^2) ./ ((v1.^2)./(n1^2*(n1-1)) + (v2.^2)./(n2^2*(n2-1)));
else
    sp2 = ((n1-1).*v1 + (n2-1).*v2) ./ (n1+n2-2);
    se2 = sp2 .* (1/n1 + 1/n2);
    t   = (m1 - m2) ./ sqrt(se2);
    df  = (n1+n2-2) * ones(size(t));
end

p = 2 * tcdf(-abs(t), df);
p(~isfinite(p)) = 1;
p = max(min(p,1),0);
end

function p = eegzoom_paired_bin_pvals(dev, stdv)
dev = double(dev);
stdv = double(stdv);

if isvector(dev), dev = dev(:); end
if isvector(stdv), stdv = stdv(:); end

n = size(dev,1);
if n < 2
    p = ones(size(dev,2),1);
    return;
end

D = dev - stdv;
m = mean(D,1,'omitnan');
s = std(D,0,1,'omitnan');

t  = m ./ (s ./ sqrt(n));
df = n - 1;

p = 2 * tcdf(-abs(t), df);
p(~isfinite(p)) = 1;
p = max(min(p,1),0);
p = p(:);
end

function sig = eegzoom_fwe_mask(pRow, alpha, method)
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

function mask = eegzoom_keep_min_run_length(mask, minLen)
mask = logical(mask(:)');
if isempty(mask) || minLen <= 1, return; end

d = diff([false mask false]);
rs = find(d==1);
re = find(d==-1)-1;

out = false(size(mask));
for i = 1:numel(rs)
    if (re(i)-rs(i)+1) >= minLen
        out(rs(i):re(i)) = true;
    end
end
mask = out;
end


function plot_topo_snaps_stats_fullscreen_nobg(EEGtmpl, tms, pCT, maskTime, cfg, bigTitle, outPng)
% Full-screen stats topoplots using the SAME p-value color style as
% your normal EEG stats topo plots.
%
% This function only plots existing pCT/maskTime.
% It does NOT change the statistics.

snap   = cfg.snapTimes_ms(:)';
snap   = snap(snap >= min(tms) & snap <= max(tms));
chlocs = EEGtmpl.chanlocs;
nSnap  = numel(snap);

if nSnap == 0
    warning('No snap times fall inside the provided time vector. No topo figure exported.');
    return;
end

% Use the same red-dot threshold logic as your normal stats topo function
dotAlpha = cfg.alpha_grp;

% Same map limits as your normal stats topo function
if isfield(cfg, 'statsTopoMapLimits') && ~isempty(cfg.statsTopoMapLimits)
    mapLimits = cfg.statsTopoMapLimits;
else
    mapLimits = [0 1];
end

% ===== SAME blue-white p-value colormap as plot_topo_snaps_stats =====
a = [247,251,255]/255;
b = [222,235,247]/255;
c = [198,219,239]/255;
d = [158,202,225]/255;
e = [107,174,214]/255;
f = [66,146,198]/255;
g = [33,113,181]/255;

customColors   = [g; f; e; d; c; b; a];
smoothGradient = interp1(linspace(0,1,size(customColors,1)), ...
    customColors, linspace(0,1,256), 'pchip');

% ===== Full-screen figure, no background =====
fig = figure( ...
    'Color', 'none', ...
    'Units', 'normalized', ...
    'OuterPosition', [0 0 1 1], ...
    'InvertHardcopy', 'off');

% Important: force figure colormap BEFORE topoplot starts
colormap(fig, smoothGradient);
set(fig, 'Colormap', smoothGradient);

% Layout: 5 + 5 + 3 maps
nCols = 5;
nRows = ceil(nSnap / nCols);

leftMargin   = 0.015;
rightMargin  = 0.015;
topMargin    = 0.085;
bottomMargin = 0.105;
gapX         = 0.015;
gapY         = 0.035;

usableW = 1 - leftMargin - rightMargin - (nCols-1)*gapX;
usableH = 1 - topMargin - bottomMargin - (nRows-1)*gapY;

axW = usableW / nCols;
axH = usableH / nRows;

for i = 1:nSnap

    row = ceil(i / nCols);
    col = i - (row-1)*nCols;

    x0 = leftMargin + (col-1) * (axW + gapX);
    y0 = 1 - topMargin - row*axH - (row-1)*gapY;

    ax = axes( ...
        'Parent', fig, ...
        'Units', 'normalized', ...
        'Position', [x0 y0 axW axH], ...
        'Color', 'none');

    % Force colormap BEFORE topoplot
    colormap(fig, smoothGradient);
    colormap(ax, smoothGradient);
    set(fig, 'Colormap', smoothGradient);

    [~, ti] = min(abs(tms - snap(i)));

    pNow = double(pCT(:,ti));
    pNow(~isfinite(pNow)) = 1;
    pNow = min(max(pNow, mapLimits(1)), mapLimits(2));

    sigCh = [];
    if logical(maskTime(ti))
        sigCh = find(pCT(:,ti) <= dotAlpha);
    end

    if isempty(sigCh)
        topoplot(pNow, chlocs, ...
            'electrodes', 'off', ...
            'maplimits', mapLimits, ...
            'numcontour', cfg.topoContours, ...
            'gridscale', cfg.topoGrid, ...
            'shading', cfg.topoShading_stats, ...
            'whitebk', 'off');
    else
        topoplot(pNow, chlocs, ...
            'electrodes', 'off', ...
            'maplimits', mapLimits, ...
            'numcontour', cfg.topoContours, ...
            'gridscale', cfg.topoGrid, ...
            'shading', cfg.topoShading_stats, ...
            'whitebk', 'off', ...
            'emarker2', {sigCh, 'o', 'r', 7, 1.2});
    end

    % CRITICAL FIX:
    % topoplot can reset the colormap internally, so force it again AFTER topoplot.
    colormap(fig, smoothGradient);
    colormap(ax, smoothGradient);
    set(fig, 'Colormap', smoothGradient);
    caxis(ax, mapLimits);

    title(ax, sprintf('%d ms', snap(i)), ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'Color', 'k');
end

% Force once more after all topoplots
colormap(fig, smoothGradient);
set(fig, 'Colormap', smoothGradient);

sgtitle(fig, bigTitle, ...
    'FontWeight', 'bold', ...
    'FontSize', 14, ...
    'Color', 'k');

% ===== Colorbar, same p-value scale =====
axCB = axes( ...
    'Parent', fig, ...
    'Units', 'normalized', ...
    'Position', [0 0 1 1], ...
    'Visible', 'off', ...
    'HandleVisibility', 'off', ...
    'Color', 'none');

colormap(axCB, smoothGradient);
set(axCB, 'CLim', mapLimits);

cb = colorbar(axCB, 'southoutside');
cb.Position = [0.72 0.045 0.25 0.025];
cb.Label.String = 'p-Values';
cb.Label.FontSize = 12;
cb.Color = 'k';
cb.Box = 'off';

% ===== Export =====
if isfield(cfg, 'exportDPI') && ~isempty(cfg.exportDPI)
    dpi = cfg.exportDPI;
else
    dpi = 300;
end

drawnow;

try
    exportgraphics(fig, outPng, ...
        'Resolution', dpi, ...
        'BackgroundColor', 'none');
catch
    warning('Transparent export failed. Exporting with white background instead.');
    exportgraphics(fig, outPng, ...
        'Resolution', dpi, ...
        'BackgroundColor', 'white');
end

close(fig);
drawnow;
end

