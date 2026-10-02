function [outSubPng, outElecPng] = eeg_export_binned_trial_matrices(outroot, grandOutFolder, cfgBIN, excludeSubs, targets)
% eeg_export_binned_trial_matrices
% PUPIL-style 2x2 matrices for EEG BINNED stats:
%   - Heatmap shows UNCORRECTED p-values (continuous)
%   - Red circles mark bins that survive FWE correction (Holm/Bonferroni)
%
% SUBJECT MATRIX (subjects x bins):
%   - Within each subject: unpaired t-tests across TRIALS, per ELECTRODE, using BIN MEANS
%   - Per-bin electrode criterion:
%         passBin(b) = (#electrodes with p<=cfgBIN.alpha_ind) >= cfgBIN.minElec_ind
%   - Heatmap value (recommended): pK(b) = k-th smallest p across electrodes (k=minElec_ind)
%       -> This directly matches the ">=k electrodes significant" rule.
%   - Red circles: FWE correction across bins on the heatmap vector pK,
%                 AND passBin(b)==true
%
% ELECTRODE MATRIX (electrodes x bins):
%   - Across subjects: paired test per ELECTRODE per BIN on (Dev-Std)
%   - Heatmap shows uncorrected p-values
%   - Red circles: FWE correction across bins within each electrode row
%
% NOTE:
%   This function is SELF-CONTAINED. Put it in its own file:
%       eeg_export_binned_trial_matrices.m
%   and call it from your EEG_stats script (see usage snippet below).
%
% Miguel Contreras-Altamirano, 2026

arguments
    outroot
    grandOutFolder
    cfgBIN struct
    excludeSubs = []
    targets = ["Frequent Local Standard","Rare Local Standard","Frequent Local Deviant","Rare Local Deviant"]
end

% normalize types (EEGLAB prefers char for filepaths)
outroot        = char(outroot);
grandOutFolder = char(grandOutFolder);
targets        = string(targets);

if iscell(excludeSubs), excludeSubs = string(excludeSubs); end
excludeSubs = string(excludeSubs);

% ---------- defaults ----------
if ~isfield(cfgBIN,'binSizeMS'), cfgBIN.binSizeMS = 100; end
if ~isfield(cfgBIN,'binStepMS'), cfgBIN.binStepMS = cfgBIN.binSizeMS; end
if ~isfield(cfgBIN,'winMS'),     cfgBIN.winMS     = [0 700]; end

if ~isfield(cfgBIN,'doFWE'),     cfgBIN.doFWE     = true; end
if ~isfield(cfgBIN,'fweMethod'), cfgBIN.fweMethod = "holm"; end   % "holm" or "bonferroni"
if ~isfield(cfgBIN,'alphaFWE'),  cfgBIN.alphaFWE  = 0.05; end

% Color scaling (match your pupil binned style; EEG often looks best with 0..0.05)
if ~isfield(cfgBIN,'climVals') || isempty(cfgBIN.climVals)
    cfgBIN.climVals = [0 0.05];
end

% Individual electrode-count criterion
if ~isfield(cfgBIN,'alpha_ind'),   cfgBIN.alpha_ind = 0.01; end
if ~isfield(cfgBIN,'minElec_ind'), cfgBIN.minElec_ind = 5; end

% Optional: mask heatmap outside passBin (default false to match pupil-style full heatmap)
if ~isfield(cfgBIN,'maskHeatmap') || isempty(cfgBIN.maskHeatmap)
    cfgBIN.maskHeatmap = false;
end

% Optional: extent across bins (default 1 = no extent)
if ~isfield(cfgBIN,'minConsecBins_ind') || isempty(cfgBIN.minConsecBins_ind)
    cfgBIN.minConsecBins_ind = 1;
end

% Trial balancing
if ~isfield(cfgBIN,'doBalance'), cfgBIN.doBalance = true; end
if ~isfield(cfgBIN,'balanceSeed'), cfgBIN.balanceSeed = 7; end

% Unpaired test type
if ~isfield(cfgBIN,'indTtest') || isempty(cfgBIN.indTtest)
    cfgBIN.indTtest = "student"; % "student" or "welch"
end

if ~exist(grandOutFolder,'dir'), mkdir(grandOutFolder); end

% ---------- find subjects ----------
dSubs = dir(fullfile(outroot,'S*'));
dSubs = dSubs([dSubs.isdir]);
names = string({dSubs.name});
names = names(~ismember(names,[".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names),'^S\d+','once'));
subNames = string(names(keep));

% numeric sort + exclude
subNums = arrayfun(@subject_num, subNames);
[~, ord] = sort(subNums);
subNames = subNames(ord);
subNums  = subNums(ord);

excNums = [];
if ~isempty(excludeSubs)
    excNums = unique(arrayfun(@subject_num, excludeSubs));
    subNames = subNames(~ismember(subNums, excNums));
    subNums  = subNums(~ismember(subNums, excNums));
end
if isempty(subNames)
    error('No participant folders found (after exclusions) under: %s', outroot);
end

% ---------- load template for times + channels ----------
EEGtmpl = [];
for s = 1:numel(subNames)
    p = subNames(s);
    fAll = fullfile(outroot, char(p), sprintf('auditory_all_LG_%s_E1.set', char(p)));
    if exist(fAll,'file')
        EEGtmpl = pop_loadset(char(fAll));
        EEGtmpl = eeg_checkset(EEGtmpl);
        break;
    end
end
if isempty(EEGtmpl)
    error('Could not find any E1 ALL set to use as template.');
end

tmplLabels = string({EEGtmpl.chanlocs.labels});
timesAll   = double(EEGtmpl.times(:))'; % ms

% bins
binEdges   = make_bin_edges(cfgBIN.winMS, cfgBIN.binSizeMS, cfgBIN.binStepMS);
binCenters = mean(binEdges, 2)'; 
nBins      = size(binEdges,1);

SUB  = struct();
ELEC = struct();
expList = ["E1","E2"];

for eidx = 1:2
    thisExp = expList(eidx);

    subjKeep = strings(0,1);

    % subject matrices
    pLoc_unc_all = [];
    pGlo_unc_all = [];
    pLoc_fwe_all = [];
    pGlo_fwe_all = [];
    passLoc_all  = [];
    passGlo_all  = [];

    % for electrode (group) matrices
    Lstd_bin_all = [];
    Ldev_bin_all = [];
    Gstd_bin_all = [];
    Gdev_bin_all = [];

    for s = 1:numel(subNames)

        participant = subNames(s);
        fAll = fullfile(outroot, char(participant), sprintf('auditory_all_LG_%s_%s.set', char(participant), char(thisExp)));
        if ~exist(fAll,'file'), continue; end

        EEGs = pop_loadset(char(fAll));
        EEGs = eeg_checkset(EEGs);
        EEGs = reorder_to_template_labels(EEGs, tmplLabels);

        if numel(EEGs.times) ~= numel(timesAll) || any(abs(double(EEGs.times(:))' - timesAll) > 1e-6)
            continue
        end
        if ~(isfield(EEGs,'epoch') && isfield(EEGs.epoch,'eventlg_type'))
            continue
        end

        % ---- 4 LG bins per trial ----
        M = false(4, EEGs.trials);
        for ep = 1:EEGs.trials
            lablist = to_strlist(EEGs.epoch(ep).eventlg_type);
            for k = 1:4
                M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
            end
        end
        idxFLS = find(M(1,:)); idxRLS = find(M(2,:)); idxFLD = find(M(3,:)); idxRLD = find(M(4,:));
        if isempty(idxFLS)||isempty(idxRLS)||isempty(idxFLD)||isempty(idxRLD), continue; end

        % factorial contrasts (same as your EEG script)
        stdLocIdx = [idxFLS(:); idxRLS(:)];
        devLocIdx = [idxFLD(:); idxRLD(:)];
        stdGloIdx = [idxFLS(:); idxFLD(:)];
        devGloIdx = [idxRLS(:); idxRLD(:)];

        % optional balancing
        if cfgBIN.doBalance
            rng(cfgBIN.balanceSeed + 1000*eidx + s);
            nMin = min(numel(stdLocIdx), numel(devLocIdx));
            stdLocIdx = stdLocIdx(randperm(numel(stdLocIdx), nMin));
            devLocIdx = devLocIdx(randperm(numel(devLocIdx), nMin));

            nMin = min(numel(stdGloIdx), numel(devGloIdx));
            stdGloIdx = stdGloIdx(randperm(numel(stdGloIdx), nMin));
            devGloIdx = devGloIdx(randperm(numel(devGloIdx), nMin));
        end

        DATA  = double(EEGs.data);
        nChan = size(DATA,1);

        % electrode x bin p-maps (within subject; across trials)
        pLoc_chbin = nan(nChan, nBins);
        pGlo_chbin = nan(nChan, nBins);

        for b = 1:nBins
            t1 = binEdges(b,1);
            t2 = binEdges(b,2);
            idxT = (timesAll >= t1) & (timesAll < t2);
            if ~any(idxT), continue; end

            % Local
            Xstd = squeeze(mean(DATA(:,idxT,stdLocIdx), 2, 'omitnan')); % [chan x trials]
            Xdev = squeeze(mean(DATA(:,idxT,devLocIdx), 2, 'omitnan')); % [chan x trials]
            pLoc_chbin(:,b) = unpaired_p_per_channel(Xdev, Xstd, cfgBIN.indTtest);

            % Global
            Xstd = squeeze(mean(DATA(:,idxT,stdGloIdx), 2, 'omitnan'));
            Xdev = squeeze(mean(DATA(:,idxT,devGloIdx), 2, 'omitnan'));
            pGlo_chbin(:,b) = unpaired_p_per_channel(Xdev, Xstd, cfgBIN.indTtest);
        end

        % =========================================================
        % SUBJECT MATRIX SCALAR p-values (match >=minElec rule)
        % =========================================================
        k = cfgBIN.minElec_ind;                         % e.g., 5
        pLoc_k = kth_smallest_per_bin(pLoc_chbin, k);   % 1 x nBins
        pGlo_k = kth_smallest_per_bin(pGlo_chbin, k);   % 1 x nBins

        % passBin = ">=k electrodes have p<=alpha_ind"
        passLoc = isfinite(pLoc_k) & (pLoc_k <= cfgBIN.alpha_ind);
        passGlo = isfinite(pGlo_k) & (pGlo_k <= cfgBIN.alpha_ind);

        % optional extent across bins (default 1 -> no extent)
        if cfgBIN.minConsecBins_ind > 1
            passLoc = apply_extent_1d(passLoc, cfgBIN.minConsecBins_ind);
            passGlo = apply_extent_1d(passGlo, cfgBIN.minConsecBins_ind);
        end

        % heatmap values (pupil-like continuous)
        pLoc_unc = pLoc_k;
        pGlo_unc = pGlo_k;

        % optional: mask outside pass (sparser, SSI-like)
        if cfgBIN.maskHeatmap
            pLoc_unc(~passLoc) = NaN;
            pGlo_unc(~passGlo) = NaN;
        end

        % FWE across bins on the SUBJECT heatmap vector
        if cfgBIN.doFWE
            pLoc_fwe = fwe_correct_vector(pLoc_unc, cfgBIN.fweMethod);
            pGlo_fwe = fwe_correct_vector(pGlo_unc, cfgBIN.fweMethod);
        else
            pLoc_fwe = pLoc_unc;
            pGlo_fwe = pGlo_unc;
        end

        % skip if nothing
        if all(isnan(pLoc_unc)) && all(isnan(pGlo_unc))
            continue
        end

        subjKeep(end+1,1)     = participant; %#ok<AGROW>
        pLoc_unc_all(end+1,:) = pLoc_unc;    %#ok<AGROW>
        pGlo_unc_all(end+1,:) = pGlo_unc;    %#ok<AGROW>
        pLoc_fwe_all(end+1,:) = pLoc_fwe;    %#ok<AGROW>
        pGlo_fwe_all(end+1,:) = pGlo_fwe;    %#ok<AGROW>
        passLoc_all(end+1,:)  = passLoc;     %#ok<AGROW>
        passGlo_all(end+1,:)  = passGlo;     %#ok<AGROW>

        % =========================================================
        % ELECTRODE MATRIX (group): store subject binned ERPs
        % =========================================================
        ERP_FLS = mean(DATA(:,:,idxFLS), 3, 'omitnan');
        ERP_RLS = mean(DATA(:,:,idxRLS), 3, 'omitnan');
        ERP_FLD = mean(DATA(:,:,idxFLD), 3, 'omitnan');
        ERP_RLD = mean(DATA(:,:,idxRLD), 3, 'omitnan');

        lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3, 'omitnan');
        ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3, 'omitnan');
        gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3, 'omitnan');
        gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3, 'omitnan');

        lstd_bin = nan(nChan,nBins); ldev_bin = nan(nChan,nBins);
        gstd_bin = nan(nChan,nBins); gdev_bin = nan(nChan,nBins);

        for b = 1:nBins
            t1 = binEdges(b,1); t2 = binEdges(b,2);
            idxT = (timesAll >= t1) & (timesAll < t2);
            if ~any(idxT), continue; end
            lstd_bin(:,b) = mean(lstd(:,idxT), 2, 'omitnan');
            ldev_bin(:,b) = mean(ldev(:,idxT), 2, 'omitnan');
            gstd_bin(:,b) = mean(gstd(:,idxT), 2, 'omitnan');
            gdev_bin(:,b) = mean(gdev(:,idxT), 2, 'omitnan');
        end

        Lstd_bin_all(end+1,:,:) = lstd_bin; %#ok<AGROW>
        Ldev_bin_all(end+1,:,:) = ldev_bin; %#ok<AGROW>
        Gstd_bin_all(end+1,:,:) = gstd_bin; %#ok<AGROW>
        Gdev_bin_all(end+1,:,:) = gdev_bin; %#ok<AGROW>
    end

    % store subject matrix
    SUB.(thisExp).subjects = subjKeep;
    SUB.(thisExp).times    = binCenters;
    SUB.(thisExp).pLoc_unc = pLoc_unc_all;
    SUB.(thisExp).pGlo_unc = pGlo_unc_all;
    SUB.(thisExp).pLoc_fwe = pLoc_fwe_all;
    SUB.(thisExp).pGlo_fwe = pGlo_fwe_all;
    SUB.(thisExp).passLoc  = passLoc_all;
    SUB.(thisExp).passGlo  = passGlo_all;

    % build electrode (group) matrices
    nSub = size(Lstd_bin_all,1);
    if nSub >= 2
        nChan = numel(tmplLabels);
        pLoc_unc_e = nan(nChan,nBins);
        pGlo_unc_e = nan(nChan,nBins);

        for ch = 1:nChan
            for b = 1:nBins
                d = squeeze(Ldev_bin_all(:,ch,b) - Lstd_bin_all(:,ch,b));
                d = d(isfinite(d));
                if numel(d) >= 2
                    pLoc_unc_e(ch,b) = p_onesample_t(d);
                end

                d = squeeze(Gdev_bin_all(:,ch,b) - Gstd_bin_all(:,ch,b));
                d = d(isfinite(d));
                if numel(d) >= 2
                    pGlo_unc_e(ch,b) = p_onesample_t(d);
                end
            end
        end

        if cfgBIN.doFWE
            pLoc_fwe_e = nan(size(pLoc_unc_e));
            pGlo_fwe_e = nan(size(pGlo_unc_e));
            for ch = 1:nChan
                pLoc_fwe_e(ch,:) = fwe_correct_vector(pLoc_unc_e(ch,:), cfgBIN.fweMethod);
                pGlo_fwe_e(ch,:) = fwe_correct_vector(pGlo_unc_e(ch,:), cfgBIN.fweMethod);
            end
        else
            pLoc_fwe_e = pLoc_unc_e;
            pGlo_fwe_e = pGlo_unc_e;
        end

        ELEC.(thisExp).times    = binCenters;
        ELEC.(thisExp).labels   = tmplLabels;
        ELEC.(thisExp).pLoc_unc = pLoc_unc_e;
        ELEC.(thisExp).pGlo_unc = pGlo_unc_e;
        ELEC.(thisExp).pLoc_fwe = pLoc_fwe_e;
        ELEC.(thisExp).pGlo_fwe = pGlo_fwe_e;
    else
        ELEC.(thisExp).times    = binCenters;
        ELEC.(thisExp).labels   = tmplLabels;
        ELEC.(thisExp).pLoc_unc = [];
        ELEC.(thisExp).pGlo_unc = [];
        ELEC.(thisExp).pLoc_fwe = [];
        ELEC.(thisExp).pGlo_fwe = [];
    end
end

outSubPng  = fullfile(grandOutFolder, sprintf('EEG_BinnedPvals_SubjectMatrix_bin%d_step%d.png', cfgBIN.binSizeMS, cfgBIN.binStepMS));
outElecPng = fullfile(grandOutFolder, sprintf('EEG_BinnedPvals_ElectrodeMatrix_bin%d_step%d.png', cfgBIN.binSizeMS, cfgBIN.binStepMS));

plot_binned_matrix_2x2_subjects(SUB, cfgBIN, outSubPng);
plot_binned_matrix_2x2_electrodes(ELEC, cfgBIN, outElecPng);

end


%% ===================== local helpers =====================

function n = subject_num(s)
s = string(s);
tok = regexp(s, '(?<=^S)\d+', 'match', 'once');
if isempty(tok), tok = regexp(s, '\d+', 'match', 'once'); end
if isempty(tok), n = NaN; else, n = str2double(tok); end
end

function EEGout = reorder_to_template_labels(EEGin, templateLabels)
curLabels = string({EEGin.chanlocs.labels});
[~, idx] = ismember(templateLabels, curLabels);
if any(idx==0)
    missing = templateLabels(idx==0);
    error('Missing channels: %s', strjoin(cellstr(missing), ', '));
end
EEGout = EEGin;
EEGout.data     = EEGin.data(idx,:,:);
EEGout.chanlocs = EEGin.chanlocs(idx);
EEGout.nbchan   = numel(idx);
EEGout = eeg_checkset(EEGout);
end

function S = to_strlist(v)
if isstring(v), S = v;
elseif ischar(v), S = string(v);
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

function binEdges = make_bin_edges(winMS, binSizeMS, binStepMS)
tStart = winMS(1);
tEnd   = winMS(2);
starts = tStart:binStepMS:(tEnd - binSizeMS);
ends   = starts + binSizeMS;
binEdges = [starts(:) ends(:)];
end

function p = unpaired_p_per_channel(Xdev, Xstd, mode)
% Xdev/Xstd: [chan x trials]
Xdev = double(Xdev); Xstd = double(Xstd);
n1 = size(Xdev,2); n2 = size(Xstd,2);

m1 = mean(Xdev,2,'omitnan'); m2 = mean(Xstd,2,'omitnan');
v1 = var(Xdev,0,2,'omitnan'); v2 = var(Xstd,0,2,'omitnan');

if strcmpi(mode,"welch")
    se2 = v1./n1 + v2./n2;
    t   = (m1-m2) ./ sqrt(se2);
    df  = (se2.^2) ./ ( (v1.^2)./(n1^2*(n1-1)) + (v2.^2)./(n2^2*(n2-1)) );
else
    sp2 = ((n1-1).*v1 + (n2-1).*v2) ./ (n1+n2-2);
    se2 = sp2 .* (1/n1 + 1/n2);
    t   = (m1-m2) ./ sqrt(se2);
    df  = (n1+n2-2) * ones(size(t));
end

p = 2 * tcdf(-abs(t), df);
p(~isfinite(p)) = 1;
p = max(min(p,1),0);
end

function p = p_onesample_t(d)
d = double(d(:));
d = d(isfinite(d));
n = numel(d);
if n < 2, p = NaN; return; end
mu = mean(d,'omitnan');
sd = std(d,0,'omitnan');
if ~isfinite(sd) || sd==0, p = 1; return; end
t = mu / (sd / sqrt(n));
df = n - 1;
p = 2 * tcdf(-abs(t), df);
p = max(min(p,1),0);
end

function pAdj = fwe_correct_vector(p, method)
method = string(lower(method));
pAdj = p;
ok = ~isnan(p) & isfinite(p);
pv = p(ok);
if isempty(pv), return; end
m = numel(pv);

switch method
    case "bonferroni"
        pvAdj = min(pv * m, 1);
    case "holm"
        [pSorted, order] = sort(pv(:), 'ascend');
        holmFactor = (m - (1:m)' + 1);
        pHolm = pSorted .* holmFactor;
        pHolm(pHolm > 1) = 1;
        pHolm = cummax(pHolm);
        pvAdj = nan(size(pv));
        pvAdj(order) = pHolm;
    otherwise
        error('Unknown FWE method: %s (use "holm" or "bonferroni")', method);
end
pAdj(ok) = pvAdj;
end

function pK = kth_smallest_per_bin(pCT, k)
% pCT: [nElec x nBins] -> pK: [1 x nBins], k-th smallest p in each bin
pCT = double(pCT);
pCT(~isfinite(pCT)) = NaN;
nBins = size(pCT,2);
pK = nan(1,nBins);
for b = 1:nBins
    v = pCT(:,b);
    v = v(isfinite(v));
    if numel(v) < k, continue; end
    v = sort(v,'ascend');
    pK(b) = v(k);
end
end

function maskExt = apply_extent_1d(mask, minConsec)
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

function plot_binned_matrix_2x2_subjects(SUB, cfgBIN, outPng)
% same palette as your pupil function
a=[247,251,255]/255; b=[222,235,247]/255; c=[198,219,239]/255;
d=[158,202,225]/255; e=[107,174,214]/255; f=[66,146,198]/255; g=[33,113,181]/255;
customColors=[g;f;e;d;c;b;a];
smoothGradient=interp1(linspace(0,1,size(customColors,1)),customColors,linspace(0,1,256),'pchip');

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.85]);
t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
t.OuterPosition = [0.05 0.06 0.86 0.90];

climVals = cfgBIN.climVals;
doMark = isfield(cfgBIN,'doFWE') && cfgBIN.doFWE;
alphaFWE = cfgBIN.alphaFWE;

    function one(ax, Punc, Pfwe, passMask, times, ttl)
        imagesc(ax, times, 1:size(Punc,1), Punc);
        set(ax,'YDir','reverse'); % 1 at top
        ylim(ax, [0.5 size(Punc,1)+0.5]);
        yticks(ax, 1:size(Punc,1));

        xlim(ax, cfgBIN.winMS);
        xt = get(ax,'XTick');
        xt = unique([cfgBIN.winMS(1) xt cfgBIN.winMS(2)]);
        set(ax,'XTick',xt);

        xlabel(ax, "Time after 5th sound onset (ms)");
        ylabel(ax, "Participants");
        title(ax, ttl, 'FontWeight','bold');

        colormap(ax, smoothGradient);
        caxis(ax, climVals);

        if doMark && ~isempty(Pfwe)
            hold(ax, 'on');
            sig = isfinite(Pfwe) & (Pfwe <= alphaFWE);
            if ~isempty(passMask)
                sig = sig & passMask; % gate dots by >=minElec at alpha_ind (and optional extent)
            end
            for r = 1:size(sig,1)
                cc = find(sig(r,:));
                if isempty(cc), continue; end
                scatter(ax, times(cc), r*ones(size(cc)), 26, 'o', ...
                    'MarkerEdgeColor','r', 'MarkerFaceColor','r');
            end
            hold(ax, 'off');
        end
    end

ax1 = nexttile(t); one(ax1, SUB.E1.pLoc_unc, SUB.E1.pLoc_fwe, SUB.E1.passLoc, SUB.E1.times, "E1 | Local (LS vs LD)");
ax2 = nexttile(t); one(ax2, SUB.E1.pGlo_unc, SUB.E1.pGlo_fwe, SUB.E1.passGlo, SUB.E1.times, "E1 | Global (GS vs GD)");
ax3 = nexttile(t); one(ax3, SUB.E2.pLoc_unc, SUB.E2.pLoc_fwe, SUB.E2.passLoc, SUB.E2.times, "E2 | Local (LS vs LD)");
ax4 = nexttile(t); one(ax4, SUB.E2.pGlo_unc, SUB.E2.pGlo_fwe, SUB.E2.passGlo, SUB.E2.times, "E2 | Global (GS vs GD)");

cbTop = colorbar(ax2); cbBot = colorbar(ax4);
cbTop.Units='normalized'; cbBot.Units='normalized';
cbTop.Position=[0.92 0.58 0.02 0.30];
cbBot.Position=[0.92 0.12 0.02 0.30];
ylabel(cbTop,"p-values"); ylabel(cbBot,"p-values");

if doMark
    note = sprintf(" | red circles: FWE-%s <= %.3f (AND >=%d elec p<=%.3g)", ...
        string(cfgBIN.fweMethod), alphaFWE, cfgBIN.minElec_ind, cfgBIN.alpha_ind);
else
    note = "";
end
sgtitle(sprintf("EEG binned trial-level t-tests per subject | win=[%d..%d] ms | bin=%d step=%d%s", ...
    cfgBIN.winMS(1), cfgBIN.winMS(2), cfgBIN.binSizeMS, cfgBIN.binStepMS, note), 'FontWeight','bold');

exportgraphics(fig, outPng, 'Resolution', 300);
close(fig);
end

function plot_binned_matrix_2x2_electrodes(ELEC, cfgBIN, outPng)
% same palette as your pupil function
a=[247,251,255]/255; b=[222,235,247]/255; c=[198,219,239]/255;
d=[158,202,225]/255; e=[107,174,214]/255; f=[66,146,198]/255; g=[33,113,181]/255;
customColors=[g;f;e;d;c;b;a];
smoothGradient=interp1(linspace(0,1,size(customColors,1)),customColors,linspace(0,1,256),'pchip');

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.85]);
t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
t.OuterPosition = [0.05 0.06 0.86 0.90];

climVals = cfgBIN.climVals;
doMark = isfield(cfgBIN,'doFWE') && cfgBIN.doFWE;
alphaFWE = cfgBIN.alphaFWE;

    function one(ax, Punc, Pfwe, times, ttl)
        if isempty(Punc)
            axis(ax,'off'); title(ax, ttl + " (no data)", 'FontWeight','bold'); return
        end
        imagesc(ax, times, 1:size(Punc,1), Punc);
        set(ax,'YDir','reverse');
        ylim(ax, [0.5 size(Punc,1)+0.5]);

        xlim(ax, cfgBIN.winMS);
        xt = get(ax,'XTick');
        xt = unique([cfgBIN.winMS(1) xt cfgBIN.winMS(2)]);
        set(ax,'XTick',xt);

        xlabel(ax, "Time after 5th sound onset (ms)");
        ylabel(ax, "Electrodes");
        title(ax, ttl, 'FontWeight','bold');

        colormap(ax, smoothGradient);
        caxis(ax, climVals);

        if doMark && ~isempty(Pfwe)
            hold(ax,'on');
            sig = isfinite(Pfwe) & (Pfwe <= alphaFWE);
            for r = 1:size(sig,1)
                cc = find(sig(r,:));
                if isempty(cc), continue; end
                scatter(ax, times(cc), r*ones(size(cc)), 18, 'o', ...
                    'MarkerEdgeColor','r', 'MarkerFaceColor','r');
            end
            hold(ax,'off');
        end
    end

ax1 = nexttile(t); one(ax1, ELEC.E1.pLoc_unc, ELEC.E1.pLoc_fwe, ELEC.E1.times, "E1 | Local (paired across subs)");
ax2 = nexttile(t); one(ax2, ELEC.E1.pGlo_unc, ELEC.E1.pGlo_fwe, ELEC.E1.times, "E1 | Global (paired across subs)");
ax3 = nexttile(t); one(ax3, ELEC.E2.pLoc_unc, ELEC.E2.pLoc_fwe, ELEC.E2.times, "E2 | Local (paired across subs)");
ax4 = nexttile(t); one(ax4, ELEC.E2.pGlo_unc, ELEC.E2.pGlo_fwe, ELEC.E2.times, "E2 | Global (paired across subs)");

cbTop = colorbar(ax2); cbBot = colorbar(ax4);
cbTop.Units='normalized'; cbBot.Units='normalized';
cbTop.Position=[0.92 0.58 0.02 0.30];
cbBot.Position=[0.92 0.12 0.02 0.30];
ylabel(cbTop,"p-values"); ylabel(cbBot,"p-values");

if doMark
    note = sprintf(" | red circles: FWE-%s <= %.3f", string(cfgBIN.fweMethod), alphaFWE);
else
    note = "";
end
sgtitle(sprintf("EEG binned GROUP paired t-tests per electrode | win=[%d..%d] ms | bin=%d step=%d%s", ...
    cfgBIN.winMS(1), cfgBIN.winMS(2), cfgBIN.binSizeMS, cfgBIN.binStepMS, note), 'FontWeight','bold');

exportgraphics(fig, outPng, 'Resolution', 300);
close(fig);
end