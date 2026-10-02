function OUT = eeg_export_interaction_effects(EEGtmpl, FLSsub, RLSsub, FLDsub, RLDsub, cfg, thisExp, outDir)
% eeg_export_interaction_effects
%
% Keeps the paper-faithful SI-style split-by-rule exports
% AND adds one TRUE OMNIBUS Local x Global interaction.
%
% INPUTS
%   EEGtmpl : template EEGLAB dataset with chanlocs and times
%   FLSsub  : [sub x chan x time] Frequent Local Standard
%   RLSsub  : [sub x chan x time] Rare     Local Standard
%   FLDsub  : [sub x chan x time] Frequent Local Deviant
%   RLDsub  : [sub x chan x time] Rare     Local Deviant
%   cfg     : config struct
%   thisExp : 'E1' or 'E2'
%   outDir  : output folder
%
% OUTPUT
%   OUT     : struct with stats + filenames
%
% Exports:
%   PAPER/SI-style outputs (unchanged):
%     1) Local effect | rule XXXXX   : RLD - FLS
%     2) Local effect | rule XXXXY   : FLD - RLS
%     3) "Interaction" Local effect x Current Rule  : (RLD - FLS) - (FLD - RLS)
%     4) Global effect | rule XXXXX  : RLD - FLS
%     5) Global effect | rule XXXXY  : RLS - FLD
%     6) "Interaction" Global effect x Current Rule : (RLD - FLS) - (RLS - FLD)
%
%   NEW TRUE OMNIBUS output:
%     7) TRUE omnibus Local x Global interaction
%        = (RLD - RLS) - (FLD - FLS)
%        = local effect at global deviant - local effect at global standard
%
% IMPORTANT:
%   Use ONLY the TRUE OMNIBUS interaction for grey-bar masks
%   in the main ERP figures.

if ~exist(outDir,'dir')
    mkdir(outDir);
end

times = double(EEGtmpl.times(:))';
tidx_plot  = find(times >= cfg.tlim_plot(1)  & times <= cfg.tlim_plot(2));
tidx_stats = find(times >= cfg.tlim_stats(1) & times <= cfg.tlim_stats(2));

tms_plot  = times(tidx_plot);
tms_stats = times(tidx_stats);

nSub = size(FLSsub,1);
if nSub < 2
    warning('Not enough subjects for interaction analysis in %s.', thisExp);
    OUT = struct();
    return;
end

expName = lgsi_expname(thisExp);

% =====================================================================
% PAPER / SUPPORTING-INFO STYLE SPLITS (KEEP AS THEY ARE)
% =====================================================================

% Local effect split by current rule
Loc_ruleXXXXX = RLDsub - FLSsub;   % RLD - FLS   (LDGD - LSGS)
Loc_ruleXXXXY = FLDsub - RLSsub;   % FLD - RLS   (LDGS - LSGD)
LocXRule      = Loc_ruleXXXXX - Loc_ruleXXXXY;

% Global effect split by current rule
Glo_ruleXXXXX = RLDsub - FLSsub;   % RLD - FLS   (LDGD - LSGS)
Glo_ruleXXXXY = RLSsub - FLDsub;   % RLS - FLD   (LSGD - LDGS)
GloXRule      = Glo_ruleXXXXX - Glo_ruleXXXXY;

% =====================================================================
% NEW: TRUE OMNIBUS LOCAL x GLOBAL INTERACTION
% =====================================================================
% Local effect at Global Standard  = FLD - FLS
% Local effect at Global Deviant   = RLD - RLS
% TRUE omnibus interaction         = (RLD - RLS) - (FLD - FLS)

LocEff_GS = FLDsub - FLSsub;   % Local effect when Global = Standard
LocEff_GD = RLDsub - RLSsub;   % Local effect when Global = Deviant
TrueOmni  = LocEff_GD - LocEff_GS;

% =====================================================================
% TITLE STRINGS
% =====================================================================

ttl_local_ruleXXXXX  = sprintf('%s | Local effect | rule XXXXX | (RLD - FLS)', expName);
ttl_local_ruleXXXXY  = sprintf('%s | Local effect | rule XXXXY | (FLD - RLS)', expName);
ttl_local_xrule      = sprintf('%s | Interaction: Local effect x Current Rule (global effect) | ((RLD - FLS) - (FLD - RLS))', expName);

ttl_global_ruleXXXXX = sprintf('%s | Global effect | rule XXXXX | (RLD - FLS)', expName);
ttl_global_ruleXXXXY = sprintf('%s | Global effect | rule XXXXY | (RLS - FLD)', expName);
ttl_global_xrule     = sprintf('%s | Interaction: Global effect x Current Rule (local effect) | ((RLD - FLS) - (RLS - FLD))', expName);

ttl_split_local  = 'Interaction | Local effect split by current rule | [RLD - FLS] vs [FLD - RLS]';
ttl_split_global = 'Interaction | Global effect split by current rule | [RLD - FLS] vs [RLS - FLD]';

ttl_true_omni_line = sprintf('%s | TRUE omnibus Local x Global interaction | [Local effect @ Global Deviant] vs [Local effect @ Global Standard]', expName);
ttl_true_omni_topo = sprintf('%s | TRUE omnibus Local x Global interaction | ((RLD - RLS) - (FLD - FLS))', expName);

% =====================================================================
% FILE NAMES
% =====================================================================

files = struct();

files.erp_local_split  = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_ERP_LOCAL_EFFECT_SPLIT_BY_RULE.png',  thisExp));
files.erp_global_split = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_ERP_GLOBAL_EFFECT_SPLIT_BY_RULE.png', thisExp));

files.local_ruleXXXXX_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_LOCAL_EFFECT_RULE_XXXXX.png',  thisExp));
files.local_ruleXXXXX_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_LOCAL_EFFECT_RULE_XXXXX.png', thisExp));

files.local_ruleXXXXY_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_LOCAL_EFFECT_RULE_XXXXY.png',  thisExp));
files.local_ruleXXXXY_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_LOCAL_EFFECT_RULE_XXXXY.png', thisExp));

files.local_xrule_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_INTERACTION_LOCAL_EFFECT_X_CURRENT_RULE_GLOBAL_EFFECT.png',  thisExp));
files.local_xrule_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_INTERACTION_LOCAL_EFFECT_X_CURRENT_RULE_GLOBAL_EFFECT.png', thisExp));

files.global_ruleXXXXX_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_GLOBAL_EFFECT_RULE_XXXXX.png',  thisExp));
files.global_ruleXXXXX_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_GLOBAL_EFFECT_RULE_XXXXX.png', thisExp));

files.global_ruleXXXXY_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_GLOBAL_EFFECT_RULE_XXXXY.png',  thisExp));
files.global_ruleXXXXY_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_GLOBAL_EFFECT_RULE_XXXXY.png', thisExp));

files.global_xrule_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_INTERACTION_GLOBAL_EFFECT_X_CURRENT_RULE_LOCAL_EFFECT.png',  thisExp));
files.global_xrule_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_INTERACTION_GLOBAL_EFFECT_X_CURRENT_RULE_LOCAL_EFFECT.png', thisExp));

% NEW TRUE OMNIBUS files
files.true_omnibus_erpLine   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_ERP_TRUE_OMNIBUS_LOCAL_GLOBAL_INTERACTION.png', thisExp));
files.true_omnibus_erpTopo   = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_ERP_TRUE_OMNIBUS_LOCAL_GLOBAL_INTERACTION.png', thisExp));
files.true_omnibus_statsTopo = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_TOPO_STATS_TRUE_OMNIBUS_LOCAL_GLOBAL_INTERACTION.png', thisExp));

files.mat = fullfile(outDir, sprintf('Bekinschtein_STATS_%s_INTERACTION_SUPPORTIVEINFO.mat', thisExp));

% =====================================================================
% ERP SPLIT LINE PLOTS (PAPER / SI STYLE)
% =====================================================================

[pLocSplit, mLocSplit, infoLocSplit] = lgsi_group_stats(Loc_ruleXXXXX(:,:,tidx_stats), Loc_ruleXXXXY(:,:,tidx_stats), cfg, 'LOCAL_RULE_SPLIT');
[pLocSplit_plot, mLocSplit_plot] = lgsi_expand_to_plot_window(pLocSplit, mLocSplit, tidx_plot, tidx_stats);

lgsi_plot_split_erps_2x2(EEGtmpl, tms_plot, ...
    Loc_ruleXXXXX(:,:,tidx_plot), Loc_ruleXXXXY(:,:,tidx_plot), ...
    pLocSplit_plot, mLocSplit_plot, cfg, thisExp, ...
    ttl_split_local, ...
    'Local effect | rule XXXXX', ...
    'Local effect | rule XXXXY', ...
    files.erp_local_split);

[pGloSplit, mGloSplit, infoGloSplit] = lgsi_group_stats(Glo_ruleXXXXX(:,:,tidx_stats), Glo_ruleXXXXY(:,:,tidx_stats), cfg, 'GLOBAL_RULE_SPLIT');
[pGloSplit_plot, mGloSplit_plot] = lgsi_expand_to_plot_window(pGloSplit, mGloSplit, tidx_plot, tidx_stats);

lgsi_plot_split_erps_2x2(EEGtmpl, tms_plot, ...
    Glo_ruleXXXXX(:,:,tidx_plot), Glo_ruleXXXXY(:,:,tidx_plot), ...
    pGloSplit_plot, mGloSplit_plot, cfg, thisExp, ...
    ttl_split_global, ...
    'Global effect | rule XXXXX', ...
    'Global effect | rule XXXXY', ...
    files.erp_global_split);

% =====================================================================
% NEW TRUE OMNIBUS ERP SPLIT LINE PLOT
% =====================================================================

[p7, m7, info7] = lgsi_group_stats(LocEff_GD(:,:,tidx_stats), LocEff_GS(:,:,tidx_stats), cfg, 'TRUE_OMNIBUS_LOCAL_GLOBAL_INTERACTION');
[p7_plot, m7_plot] = lgsi_expand_to_plot_window(p7, m7, tidx_plot, tidx_stats);

lgsi_plot_split_erps_2x2(EEGtmpl, tms_plot, ...
    LocEff_GD(:,:,tidx_plot), LocEff_GS(:,:,tidx_plot), ...
    p7_plot, m7_plot, cfg, thisExp, ...
    ttl_true_omni_line, ...
    'Local effect | Global Deviant', ...
    'Local effect | Global Standard', ...
    files.true_omnibus_erpLine);

% =====================================================================
% 1) LOCAL EFFECT | RULE XXXXX
% =====================================================================
[p1, m1, info1] = lgsi_group_stats(RLDsub(:,:,tidx_stats), FLSsub(:,:,tidx_stats), cfg, 'LOCAL_RULE_XXXXX');

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(RLDsub(:,:,tidx_stats) - FLSsub(:,:,tidx_stats),1)), cfg, ...
    ttl_local_ruleXXXXX, ...
    files.local_ruleXXXXX_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p1, m1, cfg, ...
    sprintf('%s | STATS', ttl_local_ruleXXXXX), ...
    files.local_ruleXXXXX_statsTopo);

% =====================================================================
% 2) LOCAL EFFECT | RULE XXXXY
% =====================================================================
[p2, m2, info2] = lgsi_group_stats(FLDsub(:,:,tidx_stats), RLSsub(:,:,tidx_stats), cfg, 'LOCAL_RULE_XXXXY');

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(FLDsub(:,:,tidx_stats) - RLSsub(:,:,tidx_stats),1)), cfg, ...
    ttl_local_ruleXXXXY, ...
    files.local_ruleXXXXY_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p2, m2, cfg, ...
    sprintf('%s | STATS', ttl_local_ruleXXXXY), ...
    files.local_ruleXXXXY_statsTopo);

% =====================================================================
% 3) "INTERACTION": LOCAL EFFECT x CURRENT RULE (GLOBAL EFFECT)
% =====================================================================
[p3, m3, info3] = lgsi_group_stats(Loc_ruleXXXXX(:,:,tidx_stats), Loc_ruleXXXXY(:,:,tidx_stats), cfg, 'LOCAL_X_RULE');

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(Loc_ruleXXXXX(:,:,tidx_stats) - Loc_ruleXXXXY(:,:,tidx_stats),1)), cfg, ...
    ttl_local_xrule, ...
    files.local_xrule_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p3, m3, cfg, ...
    sprintf('%s | STATS', ttl_local_xrule), ...
    files.local_xrule_statsTopo);

% =====================================================================
% 4) GLOBAL EFFECT | RULE XXXXX
% =====================================================================
[p4, m4, info4] = lgsi_group_stats(RLDsub(:,:,tidx_stats), FLSsub(:,:,tidx_stats), cfg, 'GLOBAL_RULE_XXXXX');

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(RLDsub(:,:,tidx_stats) - FLSsub(:,:,tidx_stats),1)), cfg, ...
    ttl_global_ruleXXXXX, ...
    files.global_ruleXXXXX_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p4, m4, cfg, ...
    sprintf('%s | STATS', ttl_global_ruleXXXXX), ...
    files.global_ruleXXXXX_statsTopo);

% =====================================================================
% 5) GLOBAL EFFECT | RULE XXXXY
% =====================================================================
[p5, m5, info5] = lgsi_group_stats(RLSsub(:,:,tidx_stats), FLDsub(:,:,tidx_stats), cfg, 'GLOBAL_RULE_XXXXY');

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(RLSsub(:,:,tidx_stats) - FLDsub(:,:,tidx_stats),1)), cfg, ...
    ttl_global_ruleXXXXY, ...
    files.global_ruleXXXXY_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p5, m5, cfg, ...
    sprintf('%s | STATS', ttl_global_ruleXXXXY), ...
    files.global_ruleXXXXY_statsTopo);

% =====================================================================
% 6) "INTERACTION": GLOBAL EFFECT x CURRENT RULE (LOCAL EFFECT)
% =====================================================================
[p6, m6, info6] = lgsi_group_stats(Glo_ruleXXXXX(:,:,tidx_stats), Glo_ruleXXXXY(:,:,tidx_stats), cfg, 'GLOBAL_X_RULE');

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(Glo_ruleXXXXX(:,:,tidx_stats) - Glo_ruleXXXXY(:,:,tidx_stats),1)), cfg, ...
    ttl_global_xrule, ...
    files.global_xrule_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p6, m6, cfg, ...
    sprintf('%s | STATS', ttl_global_xrule), ...
    files.global_xrule_statsTopo);

% =====================================================================
% 7) TRUE OMNIBUS LOCAL x GLOBAL INTERACTION
% =====================================================================

lgsi_plot_topo_snaps_erp(EEGtmpl, tms_stats, ...
    squeeze(mean(TrueOmni(:,:,tidx_stats),1)), cfg, ...
    ttl_true_omni_topo, ...
    files.true_omnibus_erpTopo);

lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms_stats, p7, m7, cfg, ...
    sprintf('%s | STATS', ttl_true_omni_topo), ...
    files.true_omnibus_statsTopo);

% =====================================================================
% SAVE SUMMARY
% =====================================================================
OUT = struct();
OUT.thisExp = thisExp;
OUT.nSub    = nSub;
OUT.times   = tms_stats;
OUT.files   = files;

OUT.local_ruleXXXXX.p    = p1;
OUT.local_ruleXXXXX.mask = m1;
OUT.local_ruleXXXXX.info = info1;

OUT.local_ruleXXXXY.p    = p2;
OUT.local_ruleXXXXY.mask = m2;
OUT.local_ruleXXXXY.info = info2;

OUT.local_xrule.p    = p3;
OUT.local_xrule.mask = m3;
OUT.local_xrule.info = info3;

OUT.global_ruleXXXXX.p    = p4;
OUT.global_ruleXXXXX.mask = m4;
OUT.global_ruleXXXXX.info = info4;

OUT.global_ruleXXXXY.p    = p5;
OUT.global_ruleXXXXY.mask = m5;
OUT.global_ruleXXXXY.info = info5;

OUT.global_xrule.p    = p6;
OUT.global_xrule.mask = m6;
OUT.global_xrule.info = info6;

OUT.true_omnibus.p    = p7;
OUT.true_omnibus.mask = m7;
OUT.true_omnibus.info = info7;

% THIS is the only grey-bar source to use in main ERP figures
OUT.grey_omnibus.p    = p7;
OUT.grey_omnibus.mask = m7;

OUT.split_local.info  = infoLocSplit;
OUT.split_global.info = infoGloSplit;

save(files.mat, 'OUT');

fprintf('\n[%s] Interaction export done.\n', thisExp);
disp(files);

end

%% =====================================================================
% LOCAL HELPERS FOR eeg_export_interaction_effects
% =====================================================================

function name = lgsi_expname(thisExp)
if strcmpi(thisExp,'E1')
    name = 'Counting';
else
    name = 'Visual interference';
end
end

function [pGrp, maskGrp, info] = lgsi_group_stats(A, B, cfg, tag)
pGrp = lgsi_paired_pmap(A, B);

info = struct('mcc','none','pMC',NaN,'nHit',NaN);

if lgsi_should_apply_mcc(cfg, 'group')
    [pGrp, mccInfo] = lgsi_apply_mcc_group(A, B, pGrp, cfg, tag);
    info.mcc = mccInfo;
end

maskRaw = lgsi_triple_criterion_time_mask(pGrp, cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp);
maskGrp = maskRaw;

if isfield(cfg,'mc') && isfield(cfg.mc,'doMC') && cfg.mc.doMC
    obsHas = any(maskRaw);
    if obsHas
        [pMC, nHit] = lgsi_mc_validate_triplecriterion_group( ...
            A, B, cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp, ...
            cfg.mc.nPerm, cfg.mc.seed, lgsi_getfield_default(cfg.mc,'verbose',true), ...
            lgsi_getfield_default(cfg.mc,'progressStep',100));
    else
        pMC = 1;
        nHit = cfg.mc.nPerm;
    end

    info.pMC = pMC;
    info.nHit = nHit;

    if lgsi_getfield_default(cfg.mc,'zeroMaskIfNotValidated',true)
        if pMC >= cfg.mc.alpha
            maskGrp(:) = false;
        end
    end
end
end

function p = lgsi_paired_pmap(devSub, stdSub)
D = double(devSub) - double(stdSub);
n = size(D,1);
m = squeeze(mean(D,1));
s = squeeze(std(D,0,1));
t = m ./ (s ./ sqrt(n));
df = n - 1;

p = 2 * tcdf(-abs(t), df);
p(~isfinite(p)) = 1;
p = max(min(p,1),0);
end

function maskTime = lgsi_triple_criterion_time_mask(pCT, alpha, minElec, minConsec)
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

function tf = lgsi_should_apply_mcc(cfg, level)
tf = false;
if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'applyTo') || isempty(cfg.mcc.applyTo)
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
        tf = true;
end
end

function [pAdj, info] = lgsi_apply_mcc_group(devSub, stdSub, pAlready, cfg, tag)
info = struct();

if ~isfield(cfg,'mcc') || ~isfield(cfg.mcc,'method') || isempty(cfg.mcc.method)
    info.method = 'none';
    pAdj = pAlready;
    return;
end

method = lower(string(cfg.mcc.method));
info.method = char(method);

switch method
    case "none"
        pAdj = pAlready;

    case "fdr_bh"
        q = lgsi_getfield_default(cfg.mcc,'q',0.05);
        pAdj = lgsi_fdr_bh_map(pAlready, q);
        info.q = q;

    case "holm"
        pAdj = lgsi_holm_bonferroni_map(pAlready);

    case "bonferroni"
        pAdj = lgsi_bonferroni_map(pAlready);

    case "fwe_permutation"
        fwe = cfg.mcc.fwe;
        nPerm = lgsi_getfield_default(fwe,'nPerm',2000);
        seed  = lgsi_getfield_default(fwe,'seed',777);
        verbose = lgsi_getfield_default(fwe,'verbose',true);
        step  = lgsi_getfield_default(fwe,'progressStep',200);
        [pAdj, info2] = lgsi_fwe_maxT_paired(devSub, stdSub, nPerm, seed, verbose, step);
        info = info2;
        info.method = sprintf('fwe_permutation(maxT, nPerm=%d)', nPerm);

    otherwise
        warning('[%s] Unknown MCC method: %s. Using none.', tag, method);
        pAdj = pAlready;
        info.method = 'none';
end
end

function pAdj = lgsi_fdr_bh_map(pCT, q)
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

function pAdj = lgsi_bonferroni_map(pCT)
p = pCT(:);
m = numel(p);
pAdj = min(p*m,1);
pAdj = reshape(pAdj, size(pCT));
end

function pAdj = lgsi_holm_bonferroni_map(pCT)
p = pCT(:);
m = numel(p);
[ps, idx] = sort(p,'ascend');

adj = zeros(size(ps));
for i = 1:m
    adj(i) = (m - i + 1) * ps(i);
end
adj = cummax(adj);
adj = min(adj,1);

pAdj = ones(size(p));
pAdj(idx) = adj;
pAdj = reshape(pAdj, size(pCT));
end

function [pFWE, info] = lgsi_fwe_maxT_paired(devSub, stdSub, nPerm, seed, verbose, progressStep)
rng(seed);

D = double(devSub) - double(stdSub);
nSub = size(D,1);

m = squeeze(mean(D,1));
s = squeeze(std(D,0,1));
tObs = m ./ (s ./ sqrt(nSub));
tObs(~isfinite(tObs)) = 0;
absTobs = abs(tObs(:));

maxT = zeros(nPerm,1);

for ip = 1:nPerm
    flip = rand(nSub,1) > 0.5;
    sgn = ones(nSub,1);
    sgn(flip) = -1;

    Dp = D .* reshape(sgn,[nSub 1 1]);

    mp = squeeze(mean(Dp,1));
    sp = squeeze(std(Dp,0,1));
    tp = mp ./ (sp ./ sqrt(nSub));
    tp(~isfinite(tp)) = 0;

    maxT(ip) = max(abs(tp(:)));

    if verbose && mod(ip, progressStep) == 0
        fprintf('      maxT perm %d/%d\n', ip, nPerm);
    end
end

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
end

function [pMC, nHit] = lgsi_mc_validate_triplecriterion_group(devSub, stdSub, alpha, minElec, minConsec, nPerm, seed, verbose, progressStep)
rng(seed);

nHit = 0;
nSub = size(devSub,1);

for ip = 1:nPerm
    swap = rand(nSub,1) > 0.5;

    devP = devSub;
    stdP = stdSub;

    if any(swap)
        devP(swap,:,:) = stdSub(swap,:,:);
        stdP(swap,:,:) = devSub(swap,:,:);
    end

    pPerm = lgsi_paired_pmap(devP, stdP);
    hasCrit = lgsi_triple_criterion_exists_anywhere(pPerm, alpha, minElec, minConsec);

    if hasCrit
        nHit = nHit + 1;
    end

    if verbose && mod(ip, progressStep) == 0
        fprintf('      MC perm %d/%d | hits=%d\n', ip, nPerm, nHit);
    end
end

pMC = (1 + nHit) / (1 + nPerm);
end

function tf = lgsi_triple_criterion_exists_anywhere(pCT, alpha, minElec, minConsec)
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

function [pPlot, mPlot] = lgsi_expand_to_plot_window(pStats, mStats, tidx_plot, tidx_stats)
[tf, loc] = ismember(tidx_plot, tidx_stats);

pPlot = ones(size(pStats,1), numel(tidx_plot));
mPlot = false(1, numel(tidx_plot));

pPlot(:,tf) = pStats(:,loc(tf));
mPlot(1,tf) = mStats(loc(tf));
end

function lgsi_plot_split_erps_2x2(EEGtmpl, tms, X1, X2, pMap, maskTime, cfg, thisExp, bigTitle, leg1, leg2, outPng)
chLabels = {EEGtmpl.chanlocs.labels};
fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

barCol = [0.45 0.45 0.45];

for ci = 1:4
    chanName = cfg.erpChans{ci};
    ch = find(strcmpi(chLabels, chanName), 1);
    if isempty(ch), ch = 1; end

    nexttile; hold on; grid on;

    m1 = squeeze(mean(X1(:,ch,:),1));
    m2 = squeeze(mean(X2(:,ch,:),1));

    sem1 = lgsi_sem(X1(:,ch,:));
    sem2 = lgsi_sem(X2(:,ch,:));

    lgsi_draw_sem_fill(tms, m1, sem1, [0 0.45 0.85], 0.15);
    lgsi_draw_sem_fill(tms, m2, sem2, [0.90 0.45 0.10], 0.15);

    plot(tms, m1, 'LineWidth', 2, 'Color', [0 0.45 0.85]);
    plot(tms, m2, 'LineWidth', 2, 'Color', [0.90 0.45 0.10]);

    xline(0,'k-','LineWidth',1);
    yline(0,'k-','LineWidth',1);

    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    xlim([tms(1) tms(end)]);
    ylim(cfg.ylim);
    title(sprintf('ERP at %s', upper(chanName)), 'FontWeight','bold');

    legend({leg1, leg2}, 'Location','southwest');

    tmsR = tms(:)';
    maskR = logical(maskTime(:)');
    pch = double(pMap(ch,:)); pch = pch(:)';

    showMask = maskR & (tmsR >= cfg.showFromMs) & (pch <= cfg.alpha_grp);
    showMask = lgsi_bridge_small_gaps(showMask, cfg.barMaxGap);
    showMask = lgsi_keep_min_run_length(showMask, cfg.barMinConsec);

    if any(showMask)
        yBar = cfg.ylim(2) - cfg.barYoffset_uV;

        d = diff([false showMask false]);
        rs = find(d==1);
        re = find(d==-1)-1;

        for rr = 1:numel(rs)
            idx = rs(rr):re(rr);
            plot(tmsR(idx), yBar*ones(size(idx)), '-', 'LineWidth', 6, 'Color', barCol, 'HandleVisibility','off');
        end

        if isfield(cfg,'pText') && cfg.pText && ~isempty(rs)
            text(tmsR(rs(1)), yBar-0.35, 'Interaction', ...
                'Color', barCol, 'FontSize', 9, 'FontWeight','bold', 'VerticalAlignment','top');
        end
    end
end

sgtitle(sprintf('Event Related Potentials | %s | %s', lgsi_expname(thisExp), bigTitle), 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', cfg.exportDPI);
catch
    saveas(fig, outPng);
end
close(fig); drawnow;
end

function lgsi_plot_topo_snaps_erp(EEGtmpl, tms, diffCT, cfg, bigTitle, outPng)
snap   = cfg.snapTimes_ms(:)';
chlocs = EEGtmpl.chanlocs;
nSnap  = numel(snap);

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.38]);

leftMargin   = 0.012;
rightMargin  = 0.012;
bottomMargin = 0.20;
gapX         = 0;
usableW = 1 - leftMargin - rightMargin - (nSnap-1)*gapX;
axW = usableW / nSnap;
axH = 0.58;
y0  = bottomMargin + 0.08;

ax = gobjects(nSnap,1);

for i = 1:nSnap
    x0 = leftMargin + (i-1)*(axW + gapX);
    ax(i) = axes('Parent',fig,'Units','normalized','Position',[x0 y0 axW axH]);

    [~, ti] = min(abs(tms - snap(i)));

    topoplot(double(diffCT(:,ti)), chlocs, ...
        'electrodes','on', ...
        'numcontour', cfg.topoContours, ...
        'gridscale',  cfg.topoGrid, ...
        'shading',    cfg.topoShading, ...
        'whitebk',    'off');

    if isfield(cfg,'caxis_erpTopo_interaction') && ~isempty(cfg.caxis_erpTopo_interaction)
        caxis(ax(i), cfg.caxis_erpTopo_interaction);
    elseif isfield(cfg,'caxis_erpTopo') && ~isempty(cfg.caxis_erpTopo)
        caxis(ax(i), cfg.caxis_erpTopo);
    end

    title(ax(i), sprintf('%d ms', snap(i)), 'FontSize', 12, 'FontWeight','bold');
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

function lgsi_plot_topo_snaps_stats_sameStyle(EEGtmpl, tms, pCT, maskTime, cfg, bigTitle, outPng)
snap   = cfg.snapTimes_ms(:)';
chlocs = EEGtmpl.chanlocs;
nSnap  = numel(snap);

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

leftMargin   = 0.012;
rightMargin  = 0.012;
bottomMargin = 0.20;
gapX         = 0;
usableW = 1 - leftMargin - rightMargin - (nSnap-1)*gapX;
axW = usableW / nSnap;
axH = 0.58;
y0  = bottomMargin + 0.08;

ax = gobjects(nSnap,1);

for i = 1:nSnap
    x0 = leftMargin + (i-1)*(axW + gapX);

    ax(i) = axes('Parent',fig,'Units','normalized','Position',[x0 y0 axW axH]);

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
            'emarker2',{sigCh,'o','r',5,0.8});
    end

    title(ax(i), sprintf('%d ms', snap(i)), 'FontSize', 12, 'FontWeight','bold');
end

colormap(fig, smoothGradient);

sgtitle(fig, sprintf('%s', ...
    bigTitle), 'FontWeight','bold');

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

function sem = lgsi_sem(Xsub)
X = double(squeeze(Xsub));
n = size(X,1);
if n < 2
    sem = zeros(1,size(X,2));
    return;
end
try
    s = std(X,0,1,'omitnan');
catch
    X(~isfinite(X)) = NaN;
    s = nanstd(X,0,1);
end
sem = s ./ sqrt(n);
sem = sem(:)';
end

function lgsi_draw_sem_fill(x, m, sem, faceColor, faceAlphaVal)
x   = double(x(:)');
m   = double(m(:)');
sem = double(sem(:)');

ylow  = m - sem;
yhigh = m + sem;

fill([x fliplr(x)], [ylow fliplr(yhigh)], faceColor, ...
    'FaceAlpha', faceAlphaVal, 'EdgeColor','none', 'HandleVisibility','off');
end

function mask = lgsi_keep_min_run_length(mask, minLen)
mask = logical(mask(:)');
if minLen <= 1, return; end

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

function mask = lgsi_bridge_small_gaps(mask, maxGap)
mask = logical(mask(:)');
if maxGap <= 0, return; end

d = diff([false mask false]);
rs = find(d==1);
re = find(d==-1)-1;

x = mask;
if isempty(rs), return; end

for i = 1:numel(rs)-1
    gapStart = re(i)+1;
    gapEnd   = rs(i+1)-1;
    gapLen   = gapEnd-gapStart+1;
    if gapLen > 0 && gapLen <= maxGap
        x(gapStart:gapEnd) = true;
    end
end
mask = x;
end

function v = lgsi_getfield_default(S, field, defaultVal)
if isstruct(S) && isfield(S, field)
    v = S.(field);
else
    v = defaultVal;
end
end