function lg_export_timtopo_subject_and_grand(outpath, grandFold, thisExp, targets, skipSubs, opts, EEGga)
% lg_export_timtopo_subject_and_grand
% Creates timtopo figures (2x) for each participant and for the GRAND average:
%   FIG 1 (LOCAL):  Local Standard (top) vs Local Deviant (bottom)
%   FIG 2 (GLOBAL): Global Standard (top) vs Global Deviant (bottom)
%
% Works with your auditory_all_LG_* sets and your LG labels in EEG.epoch.eventlg_type.
%
% OUTPUTS (PNG):
%   Per subject: outpath\Sxx\LG_TIMTOPO_<Exp>_LOCAL_Sxx.png
%               outpath\Sxx\LG_TIMTOPO_<Exp>_GLOBAL_Sxx.png
%   Grand avg:   grandFold\LG_TIMTOPO_<Exp>_LOCAL_GRAND_AVG.png
%               grandFold\LG_TIMTOPO_<Exp>_GLOBAL_GRAND_AVG.png
%
% Miguel Contreras-Altamirano, 2026

arguments
    outpath (1,1) string
    grandFold (1,1) string
    thisExp (1,1) string {mustBeMember(thisExp,["E1","E2"])}
    targets (1,4) string = ["Frequent Local Standard","Rare Local Standard","Frequent Local Deviant","Rare Local Deviant"]
    skipSubs string = string.empty
    opts struct = struct()
    EEGga = [] % optional: pass EEGga from your script; if empty, we try to load it from grandFold
end

% ---------------- defaults ----------------
if ~isfield(opts,'limits') || isempty(opts.limits)
    opts.limits = [-2500 700]; % show full epoch by default (matches your LG plots)
end

if ~isfield(opts,'plottimes') || isempty(opts.plottimes)
    % default requested: 13 topos between -1800 and -1300 ms
    opts.plottimes = linspace(-1800, -1300, 13);

    % Example alternative (your RP-like times):
    % opts.plottimes = [-2500 -2000 -1500 -1000 -800 -700 -600 -500 -400 -300 -200 -100 -40];
end

if ~isfield(opts,'exportDPI') || isempty(opts.exportDPI)
    opts.exportDPI = 300;
end

if ~isfield(opts,'doSubjects') || isempty(opts.doSubjects)
    opts.doSubjects = true;
end

if ~isfield(opts,'doGrand') || isempty(opts.doGrand)
    opts.doGrand = true;
end

if ~exist(grandFold,'dir'), mkdir(grandFold); end

% normalize skip list
skipSubs = string(skipSubs);
skipSubs = unique(lower(strtrim(skipSubs)));

% ---------------- find subject folders ----------------
dSubs = dir(fullfile(char(outpath),'S*'));
dSubs = dSubs([dSubs.isdir]);
names = string({dSubs.name});
names = names(~ismember(names,[".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names),'^S\d+','once'));
subNames = names(keep);

% numeric sort
subNums = arrayfun(@subject_num, subNames);
[~, ord] = sort(subNums);
subNames = subNames(ord);

% apply skip
if ~isempty(skipSubs)
    subNames = subNames(~ismember(lower(subNames), skipSubs));
end

% ---------------- template (for channel order consistency) ----------------
EEGtmpl = [];
tmplLabels = strings(0,1);

for s = 1:numel(subNames)
    p = subNames(s);
    fAll = fullfile(char(outpath), char(p), sprintf('auditory_all_LG_%s_%s.set', p, thisExp));
    if exist(fAll,'file')
        EEGtmpl = pop_loadset(char(fAll));
        EEGtmpl = eeg_checkset(EEGtmpl);
        tmplLabels = string({EEGtmpl.chanlocs.labels});
        break
    end
end

if isempty(EEGtmpl)
    warning('[%s] No subject ALL sets found -> skipping timtopo export.', thisExp);
    return;
end

% keep plottimes within limits (avoid timtopo errors)
pt = double(opts.plottimes(:))';
pt = pt(pt >= opts.limits(1) & pt <= opts.limits(2));
if numel(pt) < 2
    pt = linspace(opts.limits(1), opts.limits(2), 13);
end
opts.plottimes = pt;

% ===================== SUBJECT FIGURES =====================
if opts.doSubjects
    for s = 1:numel(subNames)
        participant = subNames(s);

        fAll = fullfile(char(outpath), char(participant), ...
            sprintf('auditory_all_LG_%s_%s.set', participant, thisExp));
        if ~exist(fAll,'file')
            continue
        end

        EEG = pop_loadset(char(fAll));
        EEG = eeg_checkset(EEG);

        % reorder channels
        EEG = reorder_to_template_labels(EEG, tmplLabels);

        % require eventlg_type
        if ~(isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type'))
            warning('[%s %s] missing EEG.epoch.eventlg_type -> skipping.', participant, thisExp);
            continue
        end

        % get 4 condition indices
        [idxFLS, idxRLS, idxFLD, idxRLD] = get_4lg_indices(EEG, targets);
        if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
            warning('[%s %s] missing one LG bin -> skipping.', participant, thisExp);
            continue
        end

        DATA = double(EEG.data);

        ERP_FLS = mean(DATA(:,:,idxFLS), 3, 'omitnan');
        ERP_RLS = mean(DATA(:,:,idxRLS), 3, 'omitnan');
        ERP_FLD = mean(DATA(:,:,idxFLD), 3, 'omitnan');
        ERP_RLD = mean(DATA(:,:,idxRLD), 3, 'omitnan');

        % factorial recombinations
        Lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3, 'omitnan');
        Ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3, 'omitnan');

        Gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3, 'omitnan');
        Gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3, 'omitnan');

        % ---- LOCAL figure ----
        outLocal = fullfile(char(outpath), char(participant), ...
            sprintf('LG_TIMTOPO_%s_LOCAL_%s.png', thisExp, participant));
        plot_timtopo_2row(Lstd, Ldev, EEG.chanlocs, opts.limits, opts.plottimes, ...
            sprintf('%s | %s | Local Standard (top) vs Local Deviant (bottom)', participant, thisExp), ...
            'Local Standard', 'Local Deviant', outLocal, opts.exportDPI);

        % ---- GLOBAL figure ----
        outGlobal = fullfile(char(outpath), char(participant), ...
            sprintf('LG_TIMTOPO_%s_GLOBAL_%s.png', thisExp, participant));
        plot_timtopo_2row(Gstd, Gdev, EEG.chanlocs, opts.limits, opts.plottimes, ...
            sprintf('%s | %s | Global Standard (top) vs Global Deviant (bottom)', participant, thisExp), ...
            'Global Standard', 'Global Deviant', outGlobal, opts.exportDPI);
    end
end

% ===================== GRAND FIGURES =====================
if opts.doGrand
    if isempty(EEGga)
        % try loading your grand set
        fGA = fullfile(char(grandFold), sprintf('Grand_avg_LG_%s.set', thisExp));
        if exist(fGA,'file')
            EEGga = pop_loadset(char(fGA));
            EEGga = eeg_checkset(EEGga);
        else
            warning('[%s] EEGga not provided and grand set not found -> skipping GRAND timtopo.', thisExp);
            return
        end
    end

    if ~(isfield(EEGga,'epoch') && isfield(EEGga.epoch,'eventlg_type'))
        warning('[%s] EEGga missing epoch.eventlg_type -> skipping GRAND timtopo.', thisExp);
        return
    end

    % compute 4 condition means across "subject-epochs"
    [idxFLS, idxRLS, idxFLD, idxRLD] = get_4lg_indices(EEGga, targets);
    if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
        warning('[%s] GRAND missing one LG bin -> skipping.', thisExp);
        return
    end

    DATAga = double(EEGga.data);
    ERP_FLS = mean(DATAga(:,:,idxFLS), 3, 'omitnan');
    ERP_RLS = mean(DATAga(:,:,idxRLS), 3, 'omitnan');
    ERP_FLD = mean(DATAga(:,:,idxFLD), 3, 'omitnan');
    ERP_RLD = mean(DATAga(:,:,idxRLD), 3, 'omitnan');

    Lstd = mean(cat(3, ERP_FLS, ERP_RLS), 3, 'omitnan');
    Ldev = mean(cat(3, ERP_FLD, ERP_RLD), 3, 'omitnan');

    Gstd = mean(cat(3, ERP_FLS, ERP_FLD), 3, 'omitnan');
    Gdev = mean(cat(3, ERP_RLS, ERP_RLD), 3, 'omitnan');

    outLocalG  = fullfile(char(grandFold), sprintf('LG_TIMTOPO_%s_LOCAL_GRAND_AVG.png', thisExp));
    outGlobalG = fullfile(char(grandFold), sprintf('LG_TIMTOPO_%s_GLOBAL_GRAND_AVG.png', thisExp));

    plot_timtopo_2row(Lstd, Ldev, EEGga.chanlocs, opts.limits, opts.plottimes, ...
        sprintf('GRAND_AVG | %s | Local Standard (top) vs Local Deviant (bottom)', thisExp), ...
        'Local Standard', 'Local Deviant', outLocalG, opts.exportDPI);

    plot_timtopo_2row(Gstd, Gdev, EEGga.chanlocs, opts.limits, opts.plottimes, ...
        sprintf('GRAND_AVG | %s | Global Standard (top) vs Global Deviant (bottom)', thisExp), ...
        'Global Standard', 'Global Deviant', outGlobalG, opts.exportDPI);
end

end % ===== end main function =====


%% ===================== helpers =====================

function n = subject_num(s)
s = string(s);
tok = regexp(s,'(?<=^S)\d+','match','once');
if isempty(tok), tok = regexp(s,'\d+','match','once'); end
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

function [idxFLS, idxRLS, idxFLD, idxRLD] = get_4lg_indices(EEG, targets)
M = false(4, EEG.trials);
for ep = 1:EEG.trials
    lablist = to_strlist(EEG.epoch(ep).eventlg_type);
    for k = 1:4
        M(k,ep) = any(strcmpi(strtrim(lablist), targets(k)));
    end
end
idxFLS = find(M(1,:));
idxRLS = find(M(2,:));
idxFLD = find(M(3,:));
idxRLD = find(M(4,:));
end

function plot_timtopo_2row(ERP_top, ERP_bot, chanlocs, limitsMS, plottimesMS, bigTitle, topTitle, botTitle, outPng, dpi)
fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]);

subplot(2,1,1);
timtopo(ERP_top, chanlocs, 'limits', limitsMS, 'plottimes', plottimesMS, 'title', topTitle);

subplot(2,1,2);
timtopo(ERP_bot, chanlocs, 'limits', limitsMS, 'plottimes', plottimesMS, 'title', botTitle);

sgtitle(bigTitle, 'FontWeight','bold');

try
    exportgraphics(fig, outPng, 'Resolution', dpi);
catch
    saveas(fig, outPng);
end
close(fig);
end