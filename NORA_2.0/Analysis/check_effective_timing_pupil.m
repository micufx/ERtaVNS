function [timing, T_LG] = check_effective_timing_pupil(participant, T_LG, thisExp, tone, varargin)
% CHECK_EFFECTIVE_TIMING_PUPIL
% Computes effective timing from EyeLink tone-onset messages:
%   SOA: within-trial diffs of tones 1..5 (4 diffs per trial)
%   ITI: TONE1(trial k) -> TONE1(trial k+1)
%
% Corrected Bekinschtein-style timing:
%   SOA = 150 ms
%   ITI = 2500-2800 ms in 50-ms steps
%
% T5 occurs 600 ms after T1, so T5(k) -> T1(k+1) is expected to be
% 1900-2200 ms. This is retained as a diagnostic but ITI_ms stores the
% pattern-onset-to-pattern-onset interval.
%
% Also stores suggested baseline window (mapped into TONE5 time) and a
% suggested epoch window start that guarantees the chosen baseline is inside.
%
% Required:
%   participant : 'S01'
%   T_LG        : Info_EEG_LG table
%   thisExp     : 'E1' or 'E2'
%   tone        : struct with tone.t_onset_ms [nTrials x 5]
%
% Stores:
%   SOA_E1_Pupil / SOA_E2_Pupil
%   ITI_E1_Pupil / ITI_E2_Pupil
%   BaselineMapped_E1_Pupil / _E2_
%   EpochSuggest_E1_Pupil / _E2_
%   BaselineAnchor_E1_Pupil / _E2_

p = inputParser;
p.addParameter('store', 'median', @(s) any(strcmpi(string(s),["median","mean"])));
p.addParameter('baselineAnchor', 'TONE5', @(s)ischar(s)||isstring(s));
p.addParameter('baselineWin_ms', [-500 0], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('epochEnd_ms', 3000, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('safety_ms', 50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('minSOA_ms', 20, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('maxSOA_ms', 600, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('minITI_ms', 2350, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('maxITI_ms', 2950, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});

storeHow       = lower(string(p.Results.store));
baselineAnchor = upper(string(p.Results.baselineAnchor));
baselineWin    = p.Results.baselineWin_ms(:)';
epochEnd_ms    = p.Results.epochEnd_ms;
safety_ms      = p.Results.safety_ms;
vb             = p.Results.verbose;

t = tone.t_onset_ms; % [nTrials x 5]

% ---- SOA (within-trial) ----
soaAll = [];
for k = 1:size(t,1)
    if any(~isfinite(t(k,:))), continue; end
    soaAll = [soaAll; diff(t(k,:)).']; %#ok<AGROW>
end
soaAll = soaAll(isfinite(soaAll));
soaAll = soaAll(soaAll>=p.Results.minSOA_ms & soaAll<=p.Results.maxSOA_ms);

% ---- ITI = T1 current -> T1 next ----
itiAll = [];
for k = 1:size(t,1)-1
    if ~isfinite(t(k,1)) || ~isfinite(t(k+1,1)), continue; end
    itiAll(end+1,1) = t(k+1,1) - t(k,1); %#ok<AGROW>
end
itiAll = itiAll(isfinite(itiAll));
itiAll = itiAll(itiAll>=p.Results.minITI_ms & itiAll<=p.Results.maxITI_ms);

% ---- Diagnostic T5 current -> T1 next ----
t5nextAll = [];
for k = 1:size(t,1)-1
    if ~isfinite(t(k,5)) || ~isfinite(t(k+1,1)), continue; end
    v = t(k+1,1) - t(k,5);
    if v>=1750 && v<=2350
        t5nextAll(end+1,1) = v; %#ok<AGROW>
    end
end

timing = struct();
timing.SOA_ms = NaN;
timing.ITI_ms = NaN;
timing.ITI_min_ms = NaN;
timing.ITI_max_ms = NaN;
timing.T5toNextT1_ms = NaN;
timing.ITI_all_ms = itiAll;
timing.T5toNextT1_all_ms = t5nextAll;

if ~isempty(soaAll)
    timing.SOA_ms = iff(storeHow=="mean", mean(soaAll,'omitnan'), median(soaAll,'omitnan'));
end
if ~isempty(itiAll)
    timing.ITI_ms = iff(storeHow=="mean", mean(itiAll,'omitnan'), median(itiAll,'omitnan'));
    timing.ITI_min_ms = min(itiAll);
    timing.ITI_max_ms = max(itiAll);
end
if ~isempty(t5nextAll)
    timing.T5toNextT1_ms = median(t5nextAll,'omitnan');
end

% ---- Suggested baseline mapped into TONE5 coordinates ----
have15 = isfinite(t(:,1)) & isfinite(t(:,5));
shift15 = t(have15,1) - t(have15,5);

if baselineAnchor=="TONE5"
    baseMapped = baselineWin;
else
    if isempty(shift15)
        baseMapped = [NaN NaN];
    else
        sh = prctile(shift15, 1);
        baseMapped = sh + baselineWin;
    end
end

epochSuggestStart = floor(baseMapped(1) - safety_ms);
epochSuggest      = [epochSuggestStart epochEnd_ms];

timing.baselineMapped_ms = baseMapped;
timing.epochSuggest_ms   = epochSuggest;

% ---- Store into T_LG ----
T_LG = ensure_pupil_timing_columns(T_LG);

rowIdx = find(strcmpi(string(T_LG.Subject_ID), string(participant)), 1);
if isempty(rowIdx)
    warning('[PupilTiming] %s not found in T_LG.Subject_ID. Not storing.', participant);
    return;
end

expStr = upper(string(thisExp));
if expStr~="E1" && expStr~="E2", expStr="E1"; end

T_LG.("SOA_"+expStr+"_Pupil")(rowIdx) = timing.SOA_ms;
T_LG.("ITI_"+expStr+"_Pupil")(rowIdx) = timing.ITI_ms;

T_LG.("BaselineMapped_"+expStr+"_Pupil"){rowIdx} = baseMapped(:).';
T_LG.("EpochSuggest_"+expStr+"_Pupil"){rowIdx}   = epochSuggest(:).';
T_LG.("BaselineAnchor_"+expStr+"_Pupil")(rowIdx) = string(baselineAnchor);

if vb
    bm = T_LG.("BaselineMapped_"+expStr+"_Pupil"){rowIdx};
    ep = T_LG.("EpochSuggest_"+expStr+"_Pupil"){rowIdx};

    fprintf(['[PupilTiming %s %s] SOA=%.0f ms | ITI T1->T1=%.0f ms ' ...
             '[%.0f %.0f] | T5->next T1=%.0f ms | baseline(mapped)=[%g %g] ' ...
             '| epochSuggest=[%g %g]\n'], ...
        participant, expStr, timing.SOA_ms, timing.ITI_ms, ...
        timing.ITI_min_ms, timing.ITI_max_ms, timing.T5toNextT1_ms, ...
        bm(1), bm(2), ep(1), ep(2));
end

function y = iff(cond, a, b)
if cond, y=a; else, y=b; end
end

end
