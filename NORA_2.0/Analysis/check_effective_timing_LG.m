function [EEG, stats, T_LG] = check_effective_timing_LG(EEG, participant, T_LG, thisExp, varargin)
% CHECK_EFFECTIVE_TIMING_LG
% Timing QC for the corrected Local-Global acquisition.
%
% Requires label_LG_events to have been run first.
%
% Expected event sequence:
%   T1 -> T2 -> T3 -> T4 -> (FLS/FLD/RLS/RLD)
%
% No marker latency correction is performed because the corrected
% experiment sends the hardware TTL at the scheduled tone onset.
%
% ITI definition used here:
%   T1 onset of pattern k -> T1 onset of pattern k+1
%
% Corrected Bekinschtein-style timing:
%   SOA = 150 ms
%   ITI = 2500-2800 ms in 50-ms steps
%
% Because T5 occurs 600 ms after T1, the corresponding
% T5(k) -> T1(k+1) interval is 1900-2200 ms.

p = inputParser;
p.addParameter('expectedSOA_ms',150,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('expectedITI_ms',2650,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('expectedITIRange_ms',[2500 2800],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('baselineBeforeX1_ms',200,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('postTone5_ms',700,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('store','median',@(s)ischar(s)||isstring(s));
p.addParameter('verbose',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('soaRange_ms',[20 600],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('itiRange_ms',[2350 2950],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('t5ToNextT1Range_ms',[1750 2350],@(x)isnumeric(x)&&numel(x)==2);
p.parse(varargin{:});
opt = p.Results;

FS = double(EEG.srate);

expTag = upper(strtrim(char(thisExp)));
if startsWith(expTag,'E')
    expNum = str2double(expTag(2:end));
else
    expNum = str2double(expTag);
end
if ~ismember(expNum,[1 2])
    expNum = 1;
end

stats = struct();
stats.SOA_ms = NaN;
stats.ITI_ms = NaN;
stats.ITI_min_ms = NaN;
stats.ITI_max_ms = NaN;
stats.T5toNextT1_ms = NaN;
stats.suggestedBaseWin_ms = [NaN NaN];
stats.suggestedEpochWin_s = [NaN NaN];
stats.nT1=0; stats.nT2=0; stats.nT3=0; stats.nT4=0; stats.nT5=0;
stats.stepMedians_ms = [NaN NaN NaN NaN];
stats.SOA_steps_all_ms = struct('d12',[],'d23',[],'d34',[],'d45',[]);
stats.ITI_all_ms = [];
stats.T5toNextT1_all_ms = [];
stats.expectedITIRange_ms = opt.expectedITIRange_ms;

if ~isfield(EEG,'event') || isempty(EEG.event)
    warning('[TIMING %s] EEG.event is empty.', expTag);
    T_LG = store_to_table(T_LG,participant,expNum,stats);
    return;
end

ev = EEG.event;
latAll = double([ev.latency]);
[latAll,ord] = sort(latAll);
ev = ev(ord);
typ = upper(strtrim(string({ev.type})));

latT1 = latAll(typ=="T1");
latT2 = latAll(typ=="T2");
latT3 = latAll(typ=="T3");
latT4 = latAll(typ=="T4");
isT5 = ismember(typ,["FLS","FLD","RLS","RLD"]);
latT5 = latAll(isT5);

stats.nT1=numel(latT1);
stats.nT2=numel(latT2);
stats.nT3=numel(latT3);
stats.nT4=numel(latT4);
stats.nT5=numel(latT5);

if opt.verbose
    fprintf('[TIMING %s] counts: T1=%d T2=%d T3=%d T4=%d T5=%d\n', ...
        expTag,stats.nT1,stats.nT2,stats.nT3,stats.nT4,stats.nT5);
end

if min([stats.nT1 stats.nT2 stats.nT3 stats.nT4 stats.nT5]) < 5
    warning('[TIMING %s] Too few markers to estimate timing robustly.', expTag);
    T_LG = store_to_table(T_LG,participant,expNum,stats);
    return;
end

n12=min(numel(latT1),numel(latT2));
n23=min(numel(latT2),numel(latT3));
n34=min(numel(latT3),numel(latT4));
n45=min(numel(latT4),numel(latT5));

d12=1000*(latT2(1:n12)-latT1(1:n12))/FS;
d23=1000*(latT3(1:n23)-latT2(1:n23))/FS;
d34=1000*(latT4(1:n34)-latT3(1:n34))/FS;
d45=1000*(latT5(1:n45)-latT4(1:n45))/FS;

d12=d12(d12>=opt.soaRange_ms(1)&d12<=opt.soaRange_ms(2));
d23=d23(d23>=opt.soaRange_ms(1)&d23<=opt.soaRange_ms(2));
d34=d34(d34>=opt.soaRange_ms(1)&d34<=opt.soaRange_ms(2));
d45=d45(d45>=opt.soaRange_ms(1)&d45<=opt.soaRange_ms(2));

stats.SOA_steps_all_ms.d12=d12;
stats.SOA_steps_all_ms.d23=d23;
stats.SOA_steps_all_ms.d34=d34;
stats.SOA_steps_all_ms.d45=d45;

stepMeds=[median(d12,'omitnan'),median(d23,'omitnan'), ...
          median(d34,'omitnan'),median(d45,'omitnan')];
stats.stepMedians_ms=stepMeds;

allSOA=[d12(:);d23(:);d34(:);d45(:)];
if strcmpi(string(opt.store),'mean')
    stats.SOA_ms=mean(allSOA,'omitnan');
else
    stats.SOA_ms=median(allSOA,'omitnan');
end

% ---- ITI = pattern onset -> next pattern onset ----
iti = 1000*diff(latT1)/FS;
iti = iti(isfinite(iti) & iti>=opt.itiRange_ms(1) & iti<=opt.itiRange_ms(2));
stats.ITI_all_ms=iti;

if ~isempty(iti)
    if strcmpi(string(opt.store),'mean')
        stats.ITI_ms=mean(iti,'omitnan');
    else
        stats.ITI_ms=median(iti,'omitnan');
    end
    stats.ITI_min_ms=min(iti);
    stats.ITI_max_ms=max(iti);
end

% Diagnostic equivalent measured from T5 of one pattern to T1 of the next.
nCross=min(numel(latT5),numel(latT1)-1);
if nCross>0
    t5next = 1000*(latT1(2:nCross+1)-latT5(1:nCross))/FS;
    t5next = t5next(isfinite(t5next) & ...
        t5next>=opt.t5ToNextT1Range_ms(1) & t5next<=opt.t5ToNextT1Range_ms(2));
    stats.T5toNextT1_all_ms=t5next;
    if ~isempty(t5next)
        stats.T5toNextT1_ms=median(t5next,'omitnan');
    end
end

offsetX1toT5_ms=sum(stepMeds,'omitnan');
baseEnd=-offsetX1toT5_ms;
baseStart=baseEnd-opt.baselineBeforeX1_ms;
stats.suggestedBaseWin_ms=[baseStart baseEnd];
stats.suggestedEpochWin_s=[baseStart/1000 opt.postTone5_ms/1000];

if opt.verbose
    fprintf('[TIMING %s] SOA(%s)=%.2f ms | steps [%.1f %.1f %.1f %.1f]\n', ...
        expTag,opt.store,stats.SOA_ms,stepMeds(1),stepMeds(2),stepMeds(3),stepMeds(4));
    fprintf('[TIMING %s] ITI T1->next T1 (%s)=%.2f ms | observed [%.1f %.1f] | expected 2500-2800 ms\n', ...
        expTag,opt.store,stats.ITI_ms,stats.ITI_min_ms,stats.ITI_max_ms);
    fprintf('[TIMING %s] T5->next T1 median=%.2f ms | expected 1900-2200 ms\n', ...
        expTag,stats.T5toNextT1_ms);
    fprintf('[TIMING %s] BaseWin_ms=[%.1f %.1f] relative to T5\n', ...
        expTag,stats.suggestedBaseWin_ms(1),stats.suggestedBaseWin_ms(2));
end

T_LG=store_to_table(T_LG,participant,expNum,stats);
EEG.etc.LG_timing=stats;

end

% =====================================================================
function T_LG=store_to_table(T_LG,participant,expNum,stats)

if nargin<1 || isempty(T_LG) || ~istable(T_LG)
    return;
end
if ~ismember('Subject_ID',T_LG.Properties.VariableNames)
    return;
end

T_LG=ensure_col(T_LG,'SOA_E1',nan(height(T_LG),1));
T_LG=ensure_col(T_LG,'SOA_E2',nan(height(T_LG),1));
T_LG=ensure_col(T_LG,'ITI_E1',nan(height(T_LG),1));
T_LG=ensure_col(T_LG,'ITI_E2',nan(height(T_LG),1));
T_LG=ensure_col(T_LG,'BaseWin_E1',repmat({[]},height(T_LG),1));
T_LG=ensure_col(T_LG,'BaseWin_E2',repmat({[]},height(T_LG),1));
T_LG=ensure_col(T_LG,'EpochWin_E1',repmat({[]},height(T_LG),1));
T_LG=ensure_col(T_LG,'EpochWin_E2',repmat({[]},height(T_LG),1));

rowIdx=find(strcmpi(string(T_LG.Subject_ID),string(participant)),1);
if isempty(rowIdx)
    return;
end

if expNum==1
    T_LG.SOA_E1(rowIdx)=stats.SOA_ms;
    T_LG.ITI_E1(rowIdx)=stats.ITI_ms;
    T_LG.BaseWin_E1{rowIdx}=stats.suggestedBaseWin_ms;
    T_LG.EpochWin_E1{rowIdx}=stats.suggestedEpochWin_s;
else
    T_LG.SOA_E2(rowIdx)=stats.SOA_ms;
    T_LG.ITI_E2(rowIdx)=stats.ITI_ms;
    T_LG.BaseWin_E2{rowIdx}=stats.suggestedBaseWin_ms;
    T_LG.EpochWin_E2{rowIdx}=stats.suggestedEpochWin_s;
end
end

% =====================================================================
function T=ensure_col(T,name,defaultValue)
if ~ismember(name,T.Properties.VariableNames)
    T.(name)=defaultValue;
end
end
