function [EEG, report] = set_LG_epoch_labels(EEG, varargin)
% SET_LG_EPOCH_LABELS
% Assign ONE unambiguous Local-Global label to each epoched EEG trial.
%
% Why this helper exists:
%   The EEG epoch is long enough to contain auditory events from neighboring
%   patterns. EEGLAB can therefore place several eventlg_* values into a
%   single EEG.epoch entry. Downstream condition masks must refer only to
%   the T5 event that actually TIME-LOCKED that epoch (latency ~= 0 ms).
%
% This function finds that zero-latency T5 anchor and writes scalar fields:
%   EEG.epoch(ep).eventlg_type
%   EEG.epoch(ep).eventlg_sound
%   EEG.epoch(ep).eventlg_effect
%   EEG.epoch(ep).eventlg_type_short
%   EEG.epoch(ep).eventlg_stim_paired / eventlg_stim_kind
%   EEG.epoch(ep).eventlg_phase / eventlg_block / eventlg_trial
%   EEG.epoch(ep).eventlg_post_stim_buffer / eventlg_exclude_primary
%
% Accepted T5 anchor types:
%   FLS, FLD, RLS, RLD
%
% It does not alter EEG.data or event latencies.

p = inputParser;
p.addParameter('tolerance_ms', [], @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
p.addParameter('verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
tol_ms = p.Results.tolerance_ms;
verbose = p.Results.verbose;

if isempty(tol_ms)
    % Slightly more than one sample.
    tol_ms = max(1, 1.5 * 1000 / double(EEG.srate));
end

report = struct();
report.nEpochs = EEG.trials;
report.nLabeled = 0;
report.nMissingAnchor = 0;
report.nMultipleAnchor = 0;
report.counts = struct('FLS',0,'FLD',0,'RLS',0,'RLD',0);

if ~isfield(EEG,'epoch') || isempty(EEG.epoch)
    warning('[LG epoch labels] EEG.epoch is empty.');
    return;
end

for ep = 1:EEG.trials

    types = toStringList(EEG.epoch(ep).eventtype);
    lats  = toNumericList(EEG.epoch(ep).eventlatency);

    n = min(numel(types), numel(lats));
    types = types(1:n);
    lats  = lats(1:n);

    isLG = ismember(upper(strtrim(types)), ["FLS","FLD","RLS","RLD"]);
    nearZero = abs(lats) <= tol_ms;
    idx = find(isLG & nearZero);

    if isempty(idx)
        % Safe fallback: among LG events, use the one closest to zero only
        % when it is still very near the epoch anchor.
        cand = find(isLG);
        if ~isempty(cand)
            [minAbs, jj] = min(abs(lats(cand)));
            if minAbs <= max(5, 2*tol_ms)
                idx = cand(jj);
            end
        end
    end

    if isempty(idx)
        report.nMissingAnchor = report.nMissingAnchor + 1;
        EEG.epoch(ep).eventlg_type = '';
        EEG.epoch(ep).eventlg_sound = '';
        EEG.epoch(ep).eventlg_effect = '';
        EEG.epoch(ep).eventlg_type_short = '';
        EEG.epoch(ep).eventlg_stim_paired = false;
        EEG.epoch(ep).eventlg_stim_kind = 'NONE';
        EEG.epoch(ep).eventlg_phase = '';
        EEG.epoch(ep).eventlg_block = NaN;
        EEG.epoch(ep).eventlg_trial = NaN;
        EEG.epoch(ep).eventlg_post_stim_buffer = false;
        EEG.epoch(ep).eventlg_exclude_primary = false;
        continue;
    end

    if numel(idx) > 1
        report.nMultipleAnchor = report.nMultipleAnchor + 1;
        [~,jj] = min(abs(lats(idx)));
        idx = idx(jj);
    end

    short = upper(strtrim(types(idx)));

    switch short
        case "FLS"
            longType = 'Frequent Local Standard';
            sound    = 'Local Standard';
            effect   = 'Global Standard';
            report.counts.FLS = report.counts.FLS + 1;
        case "FLD"
            longType = 'Frequent Local Deviant';
            sound    = 'Local Deviant';
            effect   = 'Global Standard';
            report.counts.FLD = report.counts.FLD + 1;
        case "RLS"
            longType = 'Rare Local Standard';
            sound    = 'Local Standard';
            effect   = 'Global Deviant';
            report.counts.RLS = report.counts.RLS + 1;
        case "RLD"
            longType = 'Rare Local Deviant';
            sound    = 'Local Deviant';
            effect   = 'Global Deviant';
            report.counts.RLD = report.counts.RLD + 1;
        otherwise
            continue;
    end

    EEG.epoch(ep).eventlg_type       = longType;
    EEG.epoch(ep).eventlg_sound      = sound;
    EEG.epoch(ep).eventlg_effect     = effect;
    EEG.epoch(ep).eventlg_type_short = char(short);

    % Pull stimulation/phase metadata from the actual zero-latency anchor event.
    EEG.epoch(ep).eventlg_stim_paired = false;
    EEG.epoch(ep).eventlg_stim_kind = 'NONE';
    EEG.epoch(ep).eventlg_phase = '';
    EEG.epoch(ep).eventlg_block = NaN;
    EEG.epoch(ep).eventlg_trial = NaN;
    EEG.epoch(ep).eventlg_post_stim_buffer = false;
    EEG.epoch(ep).eventlg_exclude_primary = false;
    try
        evIdxList = toNumericList(EEG.epoch(ep).event);
        if numel(evIdxList) >= idx && isfinite(evIdxList(idx))
            anchorEv = round(evIdxList(idx));
            if anchorEv>=1 && anchorEv<=numel(EEG.event)
                if isfield(EEG.event(anchorEv),'lg_stim_paired') && ~isempty(EEG.event(anchorEv).lg_stim_paired)
                    EEG.epoch(ep).eventlg_stim_paired = logical(EEG.event(anchorEv).lg_stim_paired);
                end
                if isfield(EEG.event(anchorEv),'lg_stim_kind') && ~isempty(EEG.event(anchorEv).lg_stim_kind)
                    EEG.epoch(ep).eventlg_stim_kind = EEG.event(anchorEv).lg_stim_kind;
                end
                if isfield(EEG.event(anchorEv),'lg_phase') && ~isempty(EEG.event(anchorEv).lg_phase)
                    EEG.epoch(ep).eventlg_phase = EEG.event(anchorEv).lg_phase;
                end
                if isfield(EEG.event(anchorEv),'lg_block') && ~isempty(EEG.event(anchorEv).lg_block)
                    EEG.epoch(ep).eventlg_block = EEG.event(anchorEv).lg_block;
                end
                if isfield(EEG.event(anchorEv),'lg_trial') && ~isempty(EEG.event(anchorEv).lg_trial)
                    EEG.epoch(ep).eventlg_trial = EEG.event(anchorEv).lg_trial;
                end
                if isfield(EEG.event(anchorEv),'lg_post_stim_buffer') && ~isempty(EEG.event(anchorEv).lg_post_stim_buffer)
                    EEG.epoch(ep).eventlg_post_stim_buffer = logical(EEG.event(anchorEv).lg_post_stim_buffer);
                    EEG.epoch(ep).eventlg_exclude_primary = logical(EEG.event(anchorEv).lg_post_stim_buffer);
                end
            end
        end
    catch
    end

    report.nLabeled = report.nLabeled + 1;
end

EEG.etc.LG_epoch_label_report = report;

if verbose
    fprintf('[LG epoch labels] labeled=%d/%d | FLS=%d FLD=%d RLS=%d RLD=%d\n', ...
        report.nLabeled, report.nEpochs, ...
        report.counts.FLS, report.counts.FLD, ...
        report.counts.RLS, report.counts.RLD);

    if report.nMissingAnchor > 0
        warning('[LG epoch labels] %d epoch(s) had no unique T5 anchor near 0 ms.', ...
            report.nMissingAnchor);
    end
end

end

% =====================================================================
function S = toStringList(v)
if isstring(v)
    S = v(:).';
elseif ischar(v)
    S = string({v});
elseif iscell(v)
    S = strings(1,numel(v));
    for k=1:numel(v)
        S(k)=string(v{k});
    end
elseif isnumeric(v)
    S = string(v(:).');
else
    S = strings(1,0);
end
end

% =====================================================================
function x = toNumericList(v)
if isnumeric(v)
    x = double(v(:).');
elseif iscell(v)
    x = nan(1,numel(v));
    for k=1:numel(v)
        if isnumeric(v{k}) && isscalar(v{k})
            x(k)=double(v{k});
        else
            x(k)=str2double(string(v{k}));
        end
    end
else
    x = str2double(string(v));
    x = double(x(:).');
end
end
