function plot_bekinschtein(EEG, participant, thisExp, opts)
% Local (MMN) @ topChan vs Global (P3b) @ botChan
% Uses EEG.epoch.eventlg_sound / eventlg_effect (== after your rejections)

%% ---------- defaults (NO anonymous helpers) ----------
if nargin < 4, opts = struct; end

if ~isfield(opts,'topChan')   || isempty(opts.topChan),    opts.topChan   = 'Fz';   end
if ~isfield(opts,'botChan')   || isempty(opts.botChan),    opts.botChan   = 'Pz';   end
if ~isfield(opts,'winLocal')  || isempty(opts.winLocal),   opts.winLocal  = [130 200]; end
if ~isfield(opts,'winGlobal') || isempty(opts.winGlobal),  opts.winGlobal = [320 600]; end
if ~isfield(opts,'outPng'),                                opts.outPng    = '';        end
if ~isfield(opts,'tlim')      || isempty(opts.tlim),       opts.tlim      = [EEG.times(1) EEG.times(end)]; end
if ~isfield(opts,'ylim'),                                 opts.ylim      = [];        end

times = EEG.times;

%% ---------- masks ----------
[mLocStd, mLocDev] = masks_from_epoch_field_local(EEG);
[mGloStd, mGloDev] = masks_from_epoch_field_global(EEG);

nLocStd = nnz(mLocStd); nLocDev = nnz(mLocDev);
nGloStd = nnz(mGloStd); nGloDev = nnz(mGloDev);

%% ---------- channels ----------
% logical masks
ciTop = strcmpi({EEG.chanlocs.labels}, opts.topChan);
ciBot = strcmpi({EEG.chanlocs.labels}, opts.botChan);

% fallback search
if ~any(ciTop)
    ciTop = pick_like_local({EEG.chanlocs.labels}, {'FCz','Fz','FC1','FC2'});
end
if ~any(ciBot)
    ciBot = pick_like_local({EEG.chanlocs.labels}, {'Pz','CPz','P2','CP2'});
end
% last resort
if ~any(ciTop), ciTop(1)     = true; end
if ~any(ciBot), ciBot(end)   = true; end

% convert logical -> numeric index (for markers)
idxTop = find(ciTop,1);
idxBot = find(ciBot,1);

%% ---------- ERPs ----------
Lstd = mean(double(EEG.data(idxTop,:,mLocStd)), 3);
Ldev = mean(double(EEG.data(idxTop,:,mLocDev)), 3);
Gstd = mean(double(EEG.data(idxBot,:,mGloStd)), 3);
Gdev = mean(double(EEG.data(idxBot,:,mGloDev)), 3);

%% ---------- topography data ----------
iL = times >= opts.winLocal(1)  & times <= opts.winLocal(2);
iG = times >= opts.winGlobal(1) & times <= opts.winGlobal(2);

diffLocal  = mean(double(EEG.data(:,iL,mLocDev)), [2 3]) - ...
    mean(double(EEG.data(:,iL,mLocStd)), [2 3]);
diffGlobal = mean(double(EEG.data(:,iG,mGloDev)), [2 3]) - ...
    mean(double(EEG.data(:,iG,mGloStd)), [2 3]);


% After diffLocal / diffGlobal are computed --> Color scale for topoplot
if isfield(opts,'caxis') && ~isempty(opts.caxis)
    crng = opts.caxis;                        % use caller-provided range
else
    cmax = max(abs([diffLocal; diffGlobal])); % default: from the two windows
    if cmax==0, cmax=1; end
    crng = [-cmax cmax];
end


%% ---------- figure ----------
fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w');
t   = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

%% ====== 1) Local trace ======
ax1 = nexttile(t,1); hold(ax1,'on');
ax1.FontSize = 11;   % tick numbers on both axes for the top ERP
set(groot,'defaultAxesFontSize',11);  % affects all future axes

p1 = plot(times, Lstd, 'LineWidth', 2, 'Color',[0 0.45 0.85]); %blue
p2 = plot(times, Ldev, 'LineWidth', 2, 'Color', [0.90 0.45 0.10]); %orange
yline(ax1,0,'k-'); xline(ax1,0,'k-');

if ~isempty(opts.ylim), ylim(ax1, opts.ylim); end
xlim(ax1, opts.tlim);

yl = ylim(ax1);
% shaded MMN window (clipped to tlim)
winL_clip = [max(opts.winLocal(1), opts.tlim(1)) , ...
    min(opts.winLocal(2), opts.tlim(2))];
if winL_clip(1) < winL_clip(2)
    patch(ax1, [winL_clip fliplr(winL_clip)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.82 0.90 1.00], 'FaceAlpha',0.5, 'EdgeColor','none','HandleVisibility','off');
end

grid(ax1,'on');
xlabel(ax1,'Time (ms)', FontSize=12);
ylabel(ax1,'\muV', FontSize=12);

% Local title
if isfield(EEG,'etc') && isfield(EEG.etc,'grandavg') && isfield(EEG.etc.grandavg,'nSubTotal')
    title(ax1, sprintf('%s | Local (MMN expected)   [Nsub=%d | Ntrials Std=%d Dev=%d]', ...
        upper(EEG.chanlocs(idxTop).labels), ...
        EEG.etc.grandavg.nSubTotal, ...
        EEG.etc.grandavg.nTrialLocalStd, EEG.etc.grandavg.nTrialLocalDev), ...
        'Interpreter','none', FontSize=12);
else
    title(ax1, sprintf('%s | Local (MMN expected)   [Std=%d, Dev=%d]', ...
        upper(EEG.chanlocs(idxTop).labels), nLocStd, nLocDev), ...
        'Interpreter','none', FontSize=12);
end


legend(ax1, [p1 p2], {'Local Standard','Local Deviant'}, ...
    'Location','southwest','Box','off', FontSize=12);

%% ====== 2) Local topo ======
ax2 = nexttile(t,2); cla(ax2); axes(ax2);
idxTop = find(ciTop,1);                        % numeric index

% optional: silence position warning from topoplot in tiledlayout
ws = warning('off','MATLAB:ui:Axes:UnsupportedConveniencePositionSet');

topoplot(double(diffLocal(:)), EEG.chanlocs, ...
    'electrodes','on', ...
    'emarker2',  {idxTop, 'o', [0 0 0], 12, 2}, ...
    'numcontour', 6, ...
    'gridscale', 120);

warning(ws);  % restore

% derive the shown window from the logical mask (no reliance on winLocal var)
wL = [min(times(iL)) max(times(iL))];
title(sprintf('Local at %d–%d ms', round(wL(1)), round(wL(2))), FontSize=12);
cb = colorbar(ax2,'eastoutside'); cb.Label.String = '\muV (Dev-Std)', cb.Label.FontSize = 12;

% colorbar scale topoplots:
caxis(ax2, crng);


%% ====== 3) Global trace ======
ax3 = nexttile(t,3); hold(ax3,'on');
ax3.FontSize = 11;   % bottom ERP
set(groot,'defaultAxesFontSize',11);  % affects all future axes

p3 = plot(times, Gstd, 'LineWidth', 2, 'Color',[0 0.45 0.85]);
p4 = plot(times, Gdev, 'LineWidth', 2, 'Color',[0.90 0.45 0.10]);
yline(ax3,0,'k-'); xline(ax3,0,'k-');

if ~isempty(opts.ylim), ylim(ax3, opts.ylim); end
xlim(ax3, opts.tlim);

yl = ylim(ax3);
winG_clip = [max(opts.winGlobal(1), opts.tlim(1)) , ...
    min(opts.winGlobal(2), opts.tlim(2))];
if winG_clip(1) < winG_clip(2)
    patch(ax3, [winG_clip fliplr(winG_clip)], [yl(1) yl(1) yl(2) yl(2)], ...
        [1.00 0.90 0.82], 'FaceAlpha',0.5, 'EdgeColor','none','HandleVisibility','off');
end

grid(ax3,'on');
xlabel(ax3,'Time (ms)', FontSize=12);
ylabel(ax3,'\muV', FontSize=12);



% Global title
if isfield(EEG,'etc') && isfield(EEG.etc,'grandavg') && isfield(EEG.etc.grandavg,'nSubTotal')
    title(ax3, sprintf('%s | Global (P3b expected)   [Nsub=%d | Ntrials Std=%d Dev=%d]', ...
        upper(EEG.chanlocs(idxBot).labels), ...
        EEG.etc.grandavg.nSubTotal, ...
        EEG.etc.grandavg.nTrialGlobalStd, EEG.etc.grandavg.nTrialGlobalDev), ...
        'Interpreter','none', FontSize=12);
else
    title(ax3, sprintf('%s | Global (P3b expected)   [Std=%d, Dev=%d]', ...
        upper(EEG.chanlocs(idxBot).labels), nGloStd, nGloDev), ...
        'Interpreter','none', FontSize=12);
end

legend(ax3, [p3 p4], {'Global Standard','Global Deviant'}, ...
    'Location','southwest','Box','off', FontSize=12);


%% ====== 4) Global topo ======
ax4 = nexttile(t,4); cla(ax4); axes(ax4);
idxBot = find(ciBot,1);

ws = warning('off','MATLAB:ui:Axes:UnsupportedConveniencePositionSet');

topoplot(double(diffGlobal(:)), EEG.chanlocs, ...
    'electrodes','on', ...
    'emarker2',  {idxBot, 'o', [0 0 0], 12, 2}, ...
    'numcontour', 6, ...
    'gridscale', 120);

warning(ws);

wG = [min(times(iG)) max(times(iG))];
title(sprintf('Global at %d–%d ms', round(wG(1)), round(wG(2))), FontSize=12);
cb = colorbar(ax4,'eastoutside'); cb.Label.String = '\muV (Dev-Std)', cb.Label.FontSize = 12;
% colorbar scale topoplots:
caxis(ax4, crng);


%% ====== super title ======
sgtitle(sprintf('%s | %s', participant, thisExp), ...
    'Color', "k", 'FontWeight', 'bold', 'Interpreter', 'none', FontSize=14);

%% ====== turn off toolbars (kill warning) ======
axs = findall(fig,'Type','axes');
for k = 1:numel(axs)
    try
        axtoolbar(axs(k), 'none');
    catch
    end
end


%% ====== export ======
if ~isempty(opts.outPng)
    try
        exportgraphics(fig, opts.outPng, 'Resolution', 200);
    catch
        saveas(fig, opts.outPng);
    end
end

end

%% ===== helper: local masks (safe for Local only) =====
function [mStd, mDev] = masks_from_epoch_field_local(EEG)
mStd = false(1, EEG.trials);
mDev = false(1, EEG.trials);

% 1) Preferred: explicit local labels
if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_sound')
    vals = {EEG.epoch.eventlg_sound};
    for i = 1:numel(vals)
        v = vals{i};
        % normalize to a row of strings
        if ischar(v) || isstring(v)
            vv = string(v);
        elseif iscell(v)
            vv = string(v(:)');
        else
            vv = string.empty(1,0);
        end
        mStd(i) = any(strcmpi(vv,'Local Standard'));
        mDev(i) = any(strcmpi(vv,'Local Deviant'));
    end
    return;
end

% 2) Fallback: parse Mazancieux-style combined labels in eventlg_type
if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type')
    canon = @(s) lower(strtrim(regexprep(strrep(strrep(string(s),'pattern',''),'_',' '),'\s+',' ')));
    for i = 1:numel(EEG.epoch)
        ci = "";
        try
            ci = canon(EEG.epoch(i).eventlg_type);
        catch
            % leave as "", will evaluate to false below
        end
        mStd(i) = contains(ci, "local standard");
        mDev(i) = contains(ci, "local deviant");
    end
    return;
end

% 3) If we get here, we cannot infer local labels safely
error(['masks_from_epoch_field_local: No local labels found. ' ...
       'Need EEG.epoch.eventlg_sound or eventlg_type containing "local standard/deviant".']);
end

%% ===== helper: global masks =====
function [mStd, mDev] = masks_from_epoch_field_global(EEG)
mStd = false(1, EEG.trials);
mDev = false(1, EEG.trials);

if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_effect')
    vals = {EEG.epoch.eventlg_effect};
    for i = 1:numel(vals)
        v = vals{i};
        if ischar(v) || isstring(v)
            vv = string(v);
        elseif iscell(v)
            vv = string(v(:)');
        else
            vv = string.empty(1,0);
        end
        mStd(i) = any(strcmpi(vv, 'Global Standard'));
        mDev(i) = any(strcmpi(vv, 'Global Deviant'));
    end
else
    % fallback from eventtype
    et = {EEG.epoch.eventtype};
    for i = 1:numel(et)
        v = et{i};
        if ischar(v) || isstring(v)
            vv = string(v);
        elseif iscell(v)
            vv = string(v(:)');
        else
            vv = string.empty(1,0);
        end
        mStd(i) = any(strcmpi(vv,'STANDARD'));
        mDev(i) = any(strcmpi(vv,'DEVIANT'));
    end
end
end

%% ===== helper: STANDARD/DEVIANT fallback =====
function [mStd, mDev] = masks_from_eventtype_stddev(EEG)
mStd = false(1, EEG.trials);
mDev = false(1, EEG.trials);
if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventtype')
    vals = {EEG.epoch.eventtype};
    for i = 1:numel(vals)
        v = vals{i};
        if ischar(v) || isstring(v)
            mStd(i) = strcmpi(v,'STANDARD');
            mDev(i) = strcmpi(v,'DEVIANT');
        elseif iscell(v)
            mStd(i) = any(strcmpi(v,'STANDARD'));
            mDev(i) = any(strcmpi(v,'DEVIANT'));
        end
    end
end
end

%% ===== helper: pick-like =====
function mask = pick_like_local(allLabs, prefs)
mask = false(1,numel(allLabs));
AL = upper(string(allLabs));
for k = 1:numel(prefs)
    hit = find(AL == upper(string(prefs{k})), 1);
    if ~isempty(hit)
        mask(hit) = true;
        return;
    end
end
end
