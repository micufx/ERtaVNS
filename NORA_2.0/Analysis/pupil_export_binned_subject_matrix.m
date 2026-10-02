function outBinPng = pupil_export_binned_subject_matrix(outpath, grandOutFolder, cfgBIN, excludeSubs, useField)
% pupil_export_binned_subject_matrix
% Builds per-subject binned trial-level t-tests (LS vs LD, GS vs GD) and exports a 2x2 matrix figure.
%
% If cfgBIN.doFWE == true:
%   - Heatmap shows UNCORRECTED p-values
%   - Red circles mark bins where FWE-corrected p-values <= cfgBIN.alphaFWE (default 0.05)

arguments
    outpath (1,1) string
    grandOutFolder (1,1) string
    cfgBIN struct
    excludeSubs string = []
    useField (1,1) string = "PUPIL_Z"
end

% ---- defaults ----
if ~isfield(cfgBIN,'binSizeMS'), cfgBIN.binSizeMS = 100; end
if ~isfield(cfgBIN,'binStepMS'), cfgBIN.binStepMS = cfgBIN.binSizeMS; end
if ~isfield(cfgBIN,'winMS'),     cfgBIN.winMS     = [0 2000]; end
if ~isfield(cfgBIN,'doFWE'),     cfgBIN.doFWE     = true; end
if ~isfield(cfgBIN,'fweMethod'), cfgBIN.fweMethod = "holm"; end
if ~isfield(cfgBIN,'alphaFWE'),  cfgBIN.alphaFWE  = 0.05; end

if ~exist(grandOutFolder,'dir'), mkdir(grandOutFolder); end

% ---- find subjects ----
dSubs = dir(fullfile(outpath,'S*'));
dSubs = dSubs([dSubs.isdir]);
names = string({dSubs.name});
names = names(~ismember(names,[".","..","GRAND_AVG"]));

keep = ~cellfun('isempty', regexp(cellstr(names),'^S\d+','once'));
subNames = string(names(keep));

% ---- robust numeric sorting + robust exclusion ----
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
    error('No participant folders found (after exclusions) under: %s', outpath);
end

BINNED = struct();

% ---- main loop E1/E2 ----
for eidx = 1:2
    thisExp = sprintf('E%d',eidx);

    subjKeep = strings(0,1);

    pLoc_unc_all = [];
    pGlo_unc_all = [];
    pLoc_fwe_all = [];
    pGlo_fwe_all = [];

    binCenters_ref = [];

    for s = 1:numel(subNames)

        participant = subNames(s);
        subFolder   = fullfile(outpath, participant);
        matFile     = fullfile(subFolder, sprintf('%s_%s_pupil_preproc.mat', participant, thisExp));

        if ~exist(matFile,'file'), continue; end

        S = load(matFile,'dataTableClean','relTimes');
        if ~isfield(S,'dataTableClean') || isempty(S.dataTableClean) || height(S.dataTableClean)==0
            continue
        end
        if ~isfield(S,'relTimes') || isempty(S.relTimes)
            continue
        end

        T   = S.dataTableClean;
        tms = S.relTimes(:)';

        % labels
        soundVec = string(T.sound);
        typeVec  = string(T.type);

        effectVec = strings(height(T),1);
        effectVec(startsWith(typeVec,"Frequent")) = "Global Standard";
        effectVec(startsWith(typeVec,"Rare"))     = "Global Deviant";

        % compute p-values per bin for this subject (UNC + optional FWE)
        [pLoc_unc, pLoc_fwe, binCenters] = binned_pvals_one_subject( ...
            T, useField, tms, soundVec, "Local Standard", "Local Deviant", cfgBIN);

        [pGlo_unc, pGlo_fwe, binCenters] = binned_pvals_one_subject( ...
            T, useField, tms, effectVec, "Global Standard", "Global Deviant", cfgBIN);

        if isempty(binCenters_ref)
            binCenters_ref = binCenters;
        end

        if all(isnan(pLoc_unc)) && all(isnan(pGlo_unc))
            continue
        end

        subjKeep(end+1,1) = participant; %#ok<AGROW>

        pLoc_unc_all(end+1,:) = pLoc_unc; %#ok<AGROW>
        pGlo_unc_all(end+1,:) = pGlo_unc; %#ok<AGROW>

        % keep same-size arrays always (even if doFWE=false)
        pLoc_fwe_all(end+1,:) = pLoc_fwe; %#ok<AGROW>
        pGlo_fwe_all(end+1,:) = pGlo_fwe; %#ok<AGROW>
    end

    % ---- safety: apply exclusion again at the final kept-subject list ----
    if ~isempty(excNums) && ~isempty(subjKeep)
        keepMask = ~ismember(arrayfun(@subject_num, subjKeep), excNums);
        subjKeep      = subjKeep(keepMask);
        pLoc_unc_all  = pLoc_unc_all(keepMask,:);
        pGlo_unc_all  = pGlo_unc_all(keepMask,:);
        pLoc_fwe_all  = pLoc_fwe_all(keepMask,:);
        pGlo_fwe_all  = pGlo_fwe_all(keepMask,:);
    end

    BINNED.(thisExp).subjects = subjKeep;
    BINNED.(thisExp).times    = binCenters_ref;

    BINNED.(thisExp).pLoc_unc = pLoc_unc_all;
    BINNED.(thisExp).pGlo_unc = pGlo_unc_all;

    BINNED.(thisExp).pLoc_fwe = pLoc_fwe_all;
    BINNED.(thisExp).pGlo_fwe = pGlo_fwe_all;
end

% sanity
if ~isfield(BINNED,'E1') || ~isfield(BINNED,'E2') || isempty(BINNED.E1.times)
    error('BINNED export: no data found (did not load any valid subject files).');
end

% plot + save
outBinPng = fullfile(grandOutFolder,'PUPIL_BinnedPvals_SubjectMatrix.png');
plot_binned_matrix_2x2(BINNED, cfgBIN, outBinPng);

end % ==== end MAIN function ====


%% ===================== SUBFUNCTIONS =====================

function n = subject_num(s)
s = string(s);
tok = regexp(s, '(?<=^S)\d+', 'match', 'once');
if isempty(tok)
    tok = regexp(s, '\d+', 'match', 'once');
end
if isempty(tok), n = NaN; else, n = str2double(tok); end
end


function [pUnc, pFwe, binCenters] = binned_pvals_one_subject(T, fieldName, timeVector, condVec, condA, condB, cfgBIN)
% One subject: bin trials, compute bin-mean per trial, Welch t-test between condA and condB
% Returns:
%   pUnc = uncorrected p-values per bin
%   pFwe = FWE-corrected p-values per bin (if cfgBIN.doFWE), else = pUnc

condVec = string(condVec);
idxA = (condVec == string(condA));
idxB = (condVec == string(condB));

binEdges   = make_bin_edges(cfgBIN.winMS, cfgBIN.binSizeMS, cfgBIN.binStepMS);
binCenters = mean(binEdges, 2)';         % 1 x nBins
nBins      = size(binEdges,1);

pUnc = nan(1,nBins);
pFwe = nan(1,nBins);

if ~any(idxA) || ~any(idxB)
    return
end

XA = vertcat(T.(fieldName){idxA}); % [nA x nTime]
XB = vertcat(T.(fieldName){idxB}); % [nB x nTime]

for b = 1:nBins
    t1 = binEdges(b,1);
    t2 = binEdges(b,2);

    idxT = (timeVector >= t1) & (timeVector < t2);
    if ~any(idxT)
        continue
    end

    a = mean(XA(:,idxT), 2, 'omitnan');
    c = mean(XB(:,idxT), 2, 'omitnan');

    a = a(~isnan(a));
    c = c(~isnan(c));

    if numel(a) < 2 || numel(c) < 2
        pUnc(b) = NaN;
    else
        [~, pUnc(b)] = ttest2(a, c, 'Vartype','unequal');
    end
end

% FWE correction (ONLY used for the red-circle mask)
if isfield(cfgBIN,'doFWE') && cfgBIN.doFWE
    if ~isfield(cfgBIN,'fweMethod'), cfgBIN.fweMethod = "holm"; end
    pFwe = fwe_correct_vector(pUnc, cfgBIN.fweMethod);
else
    pFwe = pUnc; % keep same size for convenience
end

end


function binEdges = make_bin_edges(winMS, binSizeMS, binStepMS)
tStart = winMS(1);
tEnd   = winMS(2);
starts = tStart:binStepMS:(tEnd - binSizeMS);
ends   = starts + binSizeMS;
binEdges = [starts(:) ends(:)];
end


function plot_binned_matrix_2x2(BINNED, cfgBIN, outPng)

% palette
a=[247,251,255]/255; b=[222,235,247]/255; c=[198,219,239]/255;
d=[158,202,225]/255; e=[107,174,214]/255; f=[66,146,198]/255; g=[33,113,181]/255;
customColors=[g;f;e;d;c;b;a];
smoothGradient=interp1(linspace(0,1,size(customColors,1)),customColors,linspace(0,1,256),'pchip');

fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.85]);

t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
t.OuterPosition = [0.05 0.06 0.86 0.90];  % leave right margin for 2 colorbars

climVals = [0 0.50]; % colormap limits for UNCORRECTED p-values

doMark = isfield(cfgBIN,'doFWE') && cfgBIN.doFWE;
if ~isfield(cfgBIN,'alphaFWE'), cfgBIN.alphaFWE = 0.05; end
alphaFWE = cfgBIN.alphaFWE;

    function one(ax, Punc, Pfwe, times, ttl)
        imagesc(ax, times, 1:size(Punc,1), Punc);
        set(ax,'YDir','reverse'); % 1 at top -> N at bottom
        ylim(ax, [0.5 size(Punc,1)+0.5]);
        yticks(ax, 1:size(Punc,1));

        % ---- NEW: force axis to include win endpoints (0 and 2000) ----
        xlim(ax, cfgBIN.winMS); % ensures endpoints are in view
        xt = get(ax,'XTick');
        xt = unique([cfgBIN.winMS(1) xt cfgBIN.winMS(2)]);
        set(ax,'XTick',xt);

        xlabel(ax, "Time after 5th sound onset (ms)");
        ylabel(ax, "Participants");
        title(ax, ttl, 'FontWeight','bold');

        colormap(ax, smoothGradient);
        caxis(ax, climVals);

        % ---- overlay red circles where FWE-corrected p <= alphaFWE ----
        if doMark && ~isempty(Pfwe)
            hold(ax, 'on');
            sig = isfinite(Pfwe) & (Pfwe <= alphaFWE);
            for r = 1:size(sig,1)
                cc = find(sig(r,:));
                if isempty(cc), continue; end
                x = times(cc);
                y = r * ones(size(x));
                scatter(ax, x, y, 26, 'o', ...
                    'MarkerEdgeColor','r', ...
                    'MarkerFaceColor','r');
            end
            hold(ax, 'off');
        end
    end

ax1 = nexttile(t); one(ax1, BINNED.E1.pLoc_unc, BINNED.E1.pLoc_fwe, BINNED.E1.times, "E1 | Local (LS vs LD)");
ax2 = nexttile(t); one(ax2, BINNED.E1.pGlo_unc, BINNED.E1.pGlo_fwe, BINNED.E1.times, "E1 | Global (GS vs GD)");
ax3 = nexttile(t); one(ax3, BINNED.E2.pLoc_unc, BINNED.E2.pLoc_fwe, BINNED.E2.times, "E2 | Local (LS vs LD)");
ax4 = nexttile(t); one(ax4, BINNED.E2.pGlo_unc, BINNED.E2.pGlo_fwe, BINNED.E2.times, "E2 | Global (GS vs GD)");

% ---- two colorbars (top row + bottom row), both for UNCORRECTED p-values ----
cbTop = colorbar(ax2);
cbBot = colorbar(ax4);

cbTop.Units = 'normalized';
cbBot.Units = 'normalized';
cbTop.Position = [0.92 0.58 0.02 0.30];
cbBot.Position = [0.92 0.12 0.02 0.30];

ylabel(cbTop, "p-values");
ylabel(cbBot, "p-values");

if doMark
    note = sprintf(" | red circles: FWE-%s <= %.3f", string(cfgBIN.fweMethod), alphaFWE);
else
    note = "";
end

sgtitle(sprintf("PUPIL binned trial-level t-tests per subject | win=[%d..%d] ms | bin=%d step=%d%s", ...
    cfgBIN.winMS(1), cfgBIN.winMS(2), cfgBIN.binSizeMS, cfgBIN.binStepMS, note), ...
    'FontWeight','bold');

exportgraphics(fig, outPng, 'Resolution', 300);
close(fig);

end


function pAdj = fwe_correct_vector(p, method)
% Adjust p-values across a vector (bins). Keeps NaNs.
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
        holmFactor = (m - (1:m)' + 1);   % m, m-1, ..., 1
        pHolm = pSorted .* holmFactor;
        pHolm(pHolm > 1) = 1;
        pHolm = cummax(pHolm);           % monotone
        pvAdj = nan(size(pv));
        pvAdj(order) = pHolm;

    otherwise
        error('Unknown FWE method: %s (use "holm" or "bonferroni")', method);
end

pAdj(ok) = pvAdj;
end