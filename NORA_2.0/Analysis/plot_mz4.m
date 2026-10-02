function plot_mz4(EEG, participant, thisExp, opts)
% Two side-by-side panels (topChan, botChan) with 4 traces (+SEM) each:
% Frequent/Rare Local Standard (blue dashed/solid) and
% Frequent/Rare Local Deviant  (orange dashed/solid).
% Labels are read from EEG.epoch.eventlg_type (post-preprocessing).

if ~exist('opts','var'), opts = struct; end
topChan = gopt(opts,'topChan','Fz');
botChan = gopt(opts,'botChan','Cz');
tlim    = gopt(opts,'tlim', [min(EEG.times) max(EEG.times)]);
ylimfix = gopt(opts,'ylim', []);     % [] -> auto
outPng  = gopt(opts,'outPng','Mazancieux');

% -------- build masks from epoch.eventlg_type (robust to multi-label epochs) --------
targets = ["Frequent Local Standard","Rare Local Standard", ...
           "Frequent Local Deviant","Rare Local Deviant"];
M = false(4, EEG.trials);

if isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type')
    for ep = 1:EEG.trials
        lablist = to_strlist(EEG.epoch(ep).eventlg_type); % -> string array of labels for this epoch
        for k = 1:4
            M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
        end
    end
else
    warning('EEG.epoch.eventlg_type missing; nothing to plot.');
end

% -------- pick channels --------
ciTop = pickChan(EEG, topChan, {'FCz','Fz','FC1','FC2'});
ciBot = pickChan(EEG, botChan, {'Pz','CPz','P2','CP2'});
labTop = EEG.chanlocs(ciTop).labels; labBot = EEG.chanlocs(ciBot).labels;

% -------- plotting style (matches Mazancieux) --------
% 1: FLS (blue dashed), 2: RLS (blue solid), 3: FLD (orange dashed), 4: RLD (orange solid)
cols   = [0 0.45 0.85; 0 0.45 0.85; 0.90 0.45 0.10; 0.90 0.45 0.10];
styles = {'--','-','--','-'};
baseLbl= {'Frequent local standard','Rare local standard', ...
          'Frequent local deviant','Rare local deviant'};

% -------- figure --------
fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w');
t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
sgtitle(sprintf('%s | %s', participant, thisExp), 'Color',"k", 'Fontweight', 'bold', 'Interpreter','none', FontSize=14);

panelPlot(EEG, ciTop, M, cols, styles, baseLbl, nexttile(t), upper(labTop), tlim, ylimfix);
panelPlot(EEG, ciBot, M, cols, styles, baseLbl, nexttile(t), upper(labBot), tlim, ylimfix);

if ~isempty(outPng)
    try, exportgraphics(fig,outPng,'Resolution',200); catch, saveas(fig,outPng); end
end
end

% ================= helpers =================
function panelPlot(EEG, ci, M, cols, styles, baseLbl, ax, chname, tlim, ylimfix)
    axes(ax); hold(ax,'on'); grid on;
    set(groot,'defaultAxesFontSize',11);
    times = EEG.times;
    set(ax,'XLim',tlim);
    if ~isempty(ylimfix), set(ax,'YLim',ylimfix); end

    legHandles = gobjects(0);

    for k = 1:4
        msk = M(k,:);
        n   = nnz(msk);
        if n==0, continue; end

        X = squeeze(double(EEG.data(ci,:,msk)));   % [1 x T x N] -> [T x N]
        if isvector(X), X = X(:); end
        if size(X,1)~=numel(times), X = X.'; end

        mu = mean(X,2);
        se = std(X,0,2)/sqrt(size(X,2));

        % SEM patch first (no legend handle)
        xx = [times, fliplr(times)];
        yy = [(mu+se).', fliplr((mu-se).')];
        fill(xx, yy, cols(k,:), 'FaceAlpha',0.15, 'EdgeColor','none', ...
            'HandleVisibility','off');

        % Line on top (goes in the legend)
        h = plot(times, mu, 'Color', cols(k,:), 'LineWidth',2, ...
                 'LineStyle', styles{k}, 'DisplayName', baseLbl{k});
        legHandles(end+1) = h; %#ok<AGROW>
    end

    xline(0,'k-'); yline(0,'k-');
    xlabel('Time (ms)', FontSize=12); 
    ylabel('\muV', FontSize=12);
    title(sprintf('%s', chname), FontSize=12);
    if ~isempty(legHandles)
        legend(legHandles, 'Location','southwest', ...
               'Box','off', 'AutoUpdate','off', FontSize=12);
    end
end

function idx = pickChan(EEG, want, fallbacks)
idx = find(strcmpi({EEG.chanlocs.labels}, want), 1);
if isempty(idx)
    for k=1:numel(fallbacks)
        idx = find(strcmpi({EEG.chanlocs.labels}, fallbacks{k}), 1);
        if ~isempty(idx), break; end
    end
end
if isempty(idx), idx = 1; end
end

function v = gopt(s, f, def)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = def; end
end

% --- robust label-to-string-array converter ---
function S = to_strlist(v)
% Returns a string array of labels for this epoch, regardless of how EEGLAB stored it
if isstring(v)
    S = v;
elseif ischar(v)
    S = string(v);
elseif iscell(v)
    % cell might contain char/strings/numerics; keep only text-like things
    S = strings(0,1);
    for i=1:numel(v)
        if isstring(v{i}) || ischar(v{i})
            S(end+1,1) = string(v{i}); %#ok<AGROW>
        end
    end
else
    S = strings(0,1); % unknown type → no labels
end
% normalize spaces a bit
S = lower(strtrim(regexprep(S,'\s+',' ')));
end
