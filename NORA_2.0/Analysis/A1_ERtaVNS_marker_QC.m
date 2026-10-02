clc; clear; close all;

%% ER-taVNS PILOT MARKER QC
% Reads a BrainVision .vmrk directly (no EEGLAB required) and verifies the
% current S1..S15 marker architecture before preprocessing.
%
% Expected auditory sequence per trial:
%   S5 -> S6 -> S7 -> S8 -> S9 -> S1..S4 -> optional S10
%
% S9 is the T5/stimulation time anchor.
% S10 only labels a stimulation-paired trial; it is NOT the physical DS5
% onset trigger.

% Timing QC settings. These are reporting thresholds, not latency corrections.
cfg.toneOnsetsMs = [-600 -450 -300 -150 0]; % T1..T5 relative to T5
cfg.stimTargetMs = 0;
cfg.toleranceMs = 5;
acquisitionMatOverride = ''; % Optional full path to the matching acquisition .mat
assert(numel(cfg.toneOnsetsMs)==5 && cfg.toneOnsetsMs(5)==0 && all(diff(cfg.toneOnsetsMs)>0),'Invalid expected tone times.');
assert(isscalar(cfg.toleranceMs) && isfinite(cfg.toleranceMs) && cfg.toleranceMs>=0,'Invalid timing tolerance.');

THIS_DIR = fileparts(mfilename('fullpath'));
PROJECT_ROOT = fileparts(THIS_DIR);
RAW_ROOT = fullfile(PROJECT_ROOT,'Experiments','Data');
OUT_ROOT = fullfile(PROJECT_ROOT,'Analysis_Output','Pilot_QC');
if exist(OUT_ROOT,'dir')~=7, mkdir(OUT_ROOT); end

answ = inputdlg({'Participant ID (e.g. S01):','Visit (1 or 2):','Experiment (E1 or E2):'}, ...
    'ER-taVNS marker QC',1,{'S01','1','E1'});
if isempty(answ), error('QC cancelled.'); end
participant = upper(strtrim(answ{1}));
visit = str2double(strtrim(answ{2}));
expName = upper(strtrim(answ{3}));
if isempty(regexp(participant,'^S\d+$','once')) || ~ismember(visit,[1 2]) || ~ismember(expName,{'E1','E2'})
    error('Use participant Sxx, visit 1/2, experiment E1/E2.');
end

stem = sprintf('%sV%d%s',participant,visit,expName);
visitFolder = fullfile(RAW_ROOT,participant,sprintf('Visit%d',visit));
vmrkPath = fullfile(visitFolder,[stem '.vmrk']);
vhdrPath = fullfile(visitFolder,[stem '.vhdr']);

if exist(vmrkPath,'file')~=2
    [f,p] = uigetfile({'*.vmrk;*.mat','BrainVision markers or acquisition MAT'},sprintf('Select BrainVision marker file for %s',stem),visitFolder);
    if isequal(f,0), error('No .vmrk selected.'); end
    [~,stem,ext] = fileparts(f);
    vmrkPath = fullfile(p,[stem '.vmrk']);
    if strcmpi(ext,'.mat'), acquisitionMatOverride=fullfile(p,f); end
    vhdrPath = fullfile(p,[stem '.vhdr']);
end

matPath=fullfile(fileparts(vmrkPath),[stem '.mat']);
if ~isempty(acquisitionMatOverride), matPath=acquisitionMatOverride; end
if ~isfile(matPath)
    [f,p]=uigetfile('*.mat','Optional: select matching acquisition MAT (Cancel = markers only)',fileparts(vmrkPath));
    if ~isequal(f,0), matPath=fullfile(p,f); end
end
QC=table(); MarkerTimingSummary=table(); fs=NaN;
if isfile(vmrkPath)
% Sampling rate from BrainVision header.
fs = NaN;
if exist(vhdrPath,'file')==2
    h = fileread(vhdrPath);
    tok = regexp(h,'(?mi)^SamplingInterval\s*=\s*([0-9.]+)','tokens','once');
    if ~isempty(tok)
        samplingInterval_us = str2double(tok{1});
        fs = 1e6/samplingInterval_us;
    end
end
if ~isfinite(fs)
    warning('Could not read SamplingInterval from .vhdr. Timing offsets will be reported in samples only.');
end

raw = fileread(vmrkPath);
lines = regexp(raw,'\r\n|\n|\r','split');
code=[]; pos=[]; mkline=[];
for i=1:numel(lines)
    s=strtrim(lines{i});
    if isempty(regexp(s,'^Mk\d+=','once')), continue; end
    rhs=regexprep(s,'^Mk\d+=','');
    parts=strsplit(rhs,',');
    if numel(parts)<3, continue; end
    desc=strtrim(parts{2});
    tok=regexp(desc,'(?i)^S\s*([0-9]+)$','tokens','once');
    if isempty(tok), continue; end
    c=str2double(tok{1});
    p0=str2double(parts{3});
    if isfinite(c)&&isfinite(p0)
        code(end+1,1)=c; %#ok<SAGROW>
        pos(end+1,1)=p0; %#ok<SAGROW>
        mkline(end+1,1)=i; %#ok<SAGROW>
    end
end

[pos,ord]=sort(pos); code=code(ord); mkline=mkline(ord);

fprintf('\n============================================================\n');
fprintf(' ER-taVNS RAW MARKER QC | %s\n',stem);
fprintf('============================================================\n');
fprintf('Marker file: %s\n',vmrkPath);
if isfinite(fs), fprintf('Sampling rate from VHDR: %.3f Hz\n',fs); end
for c=1:15
    fprintf('S%-2d count: %d\n',c,sum(code==c));
end

idxT5=find(code==9);
nT5=numel(idxT5);
trial=(1:nT5)'; block=ceil(trial/205); trialInBlock=mod(trial-1,205)+1;
phase=repmat("TEST",nT5,1); phase(trialInBlock<=25)="HABITUATION";
conditionCode=nan(nT5,1); conditionLabel=strings(nT5,1);
stimPaired=false(nT5,1); postStimBuffer=false(nT5,1); condOffsetSamples=nan(nT5,1); stimOffsetSamples=nan(nT5,1);

maxPairSamples = inf;
if isfinite(fs), maxPairSamples = round(0.050*fs); end

for t=1:nT5
    j=idxT5(t); pT5=pos(j);
    if t<nT5, pNextT5=pos(idxT5(t+1)); else, pNextT5=inf; end
    kRange=(j+1):numel(code);
    if isempty(kRange), continue; end
    kRange=kRange(pos(kRange)<pNextT5 & (pos(kRange)-pT5)<=maxPairSamples);
    kc=kRange(ismember(code(kRange),1:4));
    if ~isempty(kc)
        kc=kc(1);
        conditionCode(t)=code(kc);
        condOffsetSamples(t)=pos(kc)-pT5;
        conditionNames = {'FLS','FLD','RLS','RLD'};
        conditionLabel(t) = string(conditionNames{conditionCode(t)});
    end
    ks=kRange(code(kRange)==10);
    if ~isempty(ks)
        ks=ks(1);
        stimPaired(t)=true;
        stimOffsetSamples(t)=pos(ks)-pT5;
    end
end

isGS=ismember(conditionCode,[1 2]);
isGD=ismember(conditionCode,[3 4]);
matchedGS=stimPaired & isGS;
stimGD=stimPaired & isGD;

% Derive buffers from the deterministic sequence rule, not from a special TTL:
% every stimulation-paired trial N must be followed immediately by an
% unstimulated GS trial N+1.
ruleViolation=false(nT5,1);
for t=find(stimPaired(:))'
    if t>=nT5 || stimPaired(t+1) || ~isGS(t+1)
        ruleViolation(t)=true;
    else
        postStimBuffer(t+1)=true;
    end
end

fprintf('\nT5-linked trial decoding:\n');
fprintf('T5 anchors: %d (expected 410)\n',nT5);
fprintf('Condition labels found: %d/%d\n',sum(isfinite(conditionCode)),nT5);
fprintf('FLS=%d | FLD=%d | RLS=%d | RLD=%d (expected 169/169/36/36)\n', ...
    sum(conditionCode==1),sum(conditionCode==2),sum(conditionCode==3),sum(conditionCode==4));
fprintf('S10 paired trials=%d (expected 72)\n',sum(stimPaired));
fprintf('  matched GS=%d (expected 36)\n',sum(matchedGS));
fprintf('  stimulated GD=%d (expected 36)\n',sum(stimGD));
fprintf('  habituation paired=%d (expected 0)\n',sum(stimPaired & phase=="HABITUATION"));
fprintf('Post-stimulation GS buffers=%d (expected 72)\n',sum(postStimBuffer));
fprintf('Consecutive paired trials=%d (expected 0)\n',sum(stimPaired(1:end-1)&stimPaired(2:end)));
fprintf('N+1 buffer-rule violations=%d (expected 0)\n',sum(ruleViolation));

if isfinite(fs)
    condOffsetMs=1000*condOffsetSamples/fs;
    stimOffsetMs=1000*stimOffsetSamples/fs;
    fprintf('Median S9 -> condition offset: %.2f ms\n',median(condOffsetMs,'omitnan'));
    fprintf('Median S9 -> S10 offset:      %.2f ms\n',median(stimOffsetMs,'omitnan'));
else
    condOffsetMs=nan(size(condOffsetSamples));
    stimOffsetMs=nan(size(stimOffsetSamples));
end

QC=table(trial,block,trialInBlock,phase,pos(idxT5),conditionCode,conditionLabel, ...
    stimPaired,postStimBuffer,matchedGS,stimGD,condOffsetSamples,condOffsetMs,stimOffsetSamples,stimOffsetMs, ...
    'VariableNames',{'trial','block','trial_in_block','phase','t5_sample','condition_code', ...
    'condition','stim_paired','post_stim_buffer','matched_gs','stimulated_gd','condition_offset_samples', ...
    'condition_offset_ms','stim_offset_samples','stim_offset_ms'});

[QC,MarkerTimingSummary]=erta_tone_marker_timing(QC,code,pos,fs,cfg);
disp(MarkerTimingSummary);
outCsv=fullfile(OUT_ROOT,[stem '_marker_QC.csv']);
writetable(QC,outCsv);
fprintf('QC table saved: %s\n',outCsv);

% Hard warnings that should be resolved before using the pilot scientifically.
if nT5~=410 || sum(code==5)~=410 || sum(code==6)~=410 || sum(code==7)~=410 || sum(code==8)~=410
    warning('Auditory tone counts differ from the expected 410 per tone position.');
end
if any([sum(conditionCode==1),sum(conditionCode==2),sum(conditionCode==3),sum(conditionCode==4)] ~= [169 169 36 36])
    warning('Condition counts differ from expected FLS/FLD/RLS/RLD = 169/169/36/36.');
end
if sum(stimPaired)~=72 || sum(matchedGS)~=36 || sum(stimGD)~=36
    warning('Paired-stimulation counts differ from expected 72 total = 36 GS + 36 GD.');
end
if any(stimPaired & phase=="HABITUATION")
    warning('At least one habituation trial has S10. This should not happen.');
end
if any(stimPaired(1:end-1)&stimPaired(2:end))
    warning('Consecutive stimulation-paired trials were detected.');
end
if sum(postStimBuffer)~=72 || any(ruleViolation)
    warning('Mandatory N+1 post-stimulation GS buffer rule failed.');
end
if any(postStimBuffer & ~isGS) || any(postStimBuffer & stimPaired)
    warning('A post-stimulation buffer is not an unstimulated GS.');
end

fprintf('============================================================\n');
fprintf('Remember: S10 is a paired-trial LABEL. S9 is the T5 EEG anchor.\n');
fprintf('============================================================\n');

else
    fprintf('No marker file: acquisition-log timing only.\n');
end

TimingTrials=table(); LogTimingSummary=table(); MissingFields=table(); experimentInfo=struct();
if isfile(matPath)
    [~,matStem]=fileparts(matPath);
    assert(strcmpi(matStem,stem),'MAT and marker stems differ; select the same run.');
    [TimingTrials,LogTimingSummary,MissingFields,experimentInfo]=erta_log_timing(matPath,cfg);
    if isfield(experimentInfo,'fileStem')
        assert(strcmpi(string(experimentInfo.fileStem),stem),'MAT experimentInfo.fileStem differs from requested run.');
    end
    writetable(TimingTrials,fullfile(OUT_ROOT,[stem '_onset_timing_trials.csv']));
    writetable(MissingFields,fullfile(OUT_ROOT,[stem '_onset_missing_fields.csv']));
    % Independent reports: never silently join rows after a missing/extra EEG marker.
    if ~isempty(QC) && height(QC)~=height(TimingTrials)
        warning('Marker/log trial counts differ (%d/%d); no rowwise alignment was assumed.',height(QC),height(TimingTrials));
    end
else
    warning('Acquisition MAT missing. Marker intervals cannot establish stimulation pulse onset.');
end
TimingSummary=[MarkerTimingSummary;LogTimingSummary];
if ~isempty(TimingSummary)
    writetable(TimingSummary,fullfile(OUT_ROOT,[stem '_onset_timing_summary.csv']));
end
save(fullfile(OUT_ROOT,[stem '_onset_timing_QC.mat']),'QC','TimingTrials','TimingSummary','MissingFields','experimentInfo','cfg','matPath','vmrkPath');
reportPath=fullfile(OUT_ROOT,[stem '_onset_timing_report.txt']);
fid=fopen(reportPath,'w'); assert(fid>=0,'Cannot create timing report.');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'ER-taVNS onset timing QC: %s\nPositive offsets = later than reference. All summaries in milliseconds.\n',stem);
fprintf(fid,'MAT: %s\nMarkers: %s\n',matPath,vmrkPath);
fprintf(fid,'Software stim_trigger_getsecs = expected pulse onset, not measured current.\n');
fprintf(fid,'Actual T5 is a playback-derived timestamp, not an independent microphone measurement.\n');
fprintf(fid,'S10 in task runs labels a paired trial; its delay is NOT the stimulation delay.\n');
fprintf(fid,'T1-T4 accuracy is relative to T5 in the EEG marker clock; no absolute clock alignment is assumed.\n');
fprintf(fid,'Absent, empty, NaN and Inf timing fields remain missing; they are never replaced with zero.\n\n');
fprintf(fid,'%s\n',evalc('disp(TimingSummary)'));
if ~isempty(TimingTrials)
    fprintf(fid,'Stimulated trials: %d; valid expected offsets: %d; near -400 ms (+/-50): %d; outside tolerance: %d\n', ...
        nnz(TimingTrials.is_stimulated==1),nnz(TimingTrials.stim_timing_valid), ...
        nnz(TimingTrials.stim_near_minus400),nnz(TimingTrials.stim_outside_tolerance));
    fprintf(fid,'\nMissing fields (all and stimulated-only counts):\n%s\n',evalc('disp(MissingFields)'));
end
clear cleanup;
fprintf('Onset timing tables and report saved in: %s\n',OUT_ROOT);

function [T,S,F,info]=erta_log_timing(matPath,cfg)
% Read-only acquisition-log inspection. GetSecs and EEG sample clocks are NOT mixed.
D=load(matPath); names=fieldnames(D);
hit=names(startsWith(names,'events_audio_'));
assert(numel(hit)==1,'Expected exactly one events_audio_* variable in %s.',matPath);
E=D.(hit{1}); assert(isstruct(E) && ~isempty(E),'Empty/non-struct auditory trial log.');
info=struct(); if isfield(D,'experimentInfo'), info=D.experimentInfo; end
N=numel(E); T=table((1:N)','VariableNames',{'log_row'});
numericFields={'session_trial','block','trial_in_block','participant','visit','is_habituation', ...
    'is_stimulated','stim_hardware_available','stim_marker_code','stim_expected_burst_ms', ...
    'stim_arm_getsecs','stim_software_start_getsecs','stim_trigger_getsecs', ...
    'stim_t5_offset_ms','stim_actual_t5_offset_ms','scheduled_pattern_onset', ...
    'actual_pattern_onset','scheduled_fifth_tone_onset','actual_fifth_tone_onset','t5_ttl_write_getsecs'};
F=table();
for j=1:numel(numericFields)
    name=numericFields{j}; [v,absent,invalid]=erta_numeric_field(E,name); T.(name)=v;
    F=[F;table(string(name),sum(absent),sum(isnan(v)&~absent&~invalid),sum(isinf(v)),sum(invalid), ...
        'VariableNames',{'field','absent','NaN','Inf','invalid_or_empty'})]; %#ok<AGROW>
end
for name={'phase','physical_pattern','type_short','local_short','global_short','stim_condition','stim_site','stim_delivery_mode'}
    v=strings(N,1); v(:)=missing;
    if isfield(E,name{1})
        for j=1:N
            z=E(j).(name{1});
            if (ischar(z) && isrow(z)) || (isstring(z) && isscalar(z)), v(j)=string(z); end
        end
    end
    T.(name{1})=v;
end
T.expected_stim_vs_scheduled_t5_ms=1000*(T.stim_trigger_getsecs-T.scheduled_fifth_tone_onset);
T.expected_stim_vs_actual_t5_ms=1000*(T.stim_trigger_getsecs-T.actual_fifth_tone_onset);
T.software_start_vs_t5_ms=1000*(T.stim_software_start_getsecs-T.scheduled_fifth_tone_onset);
T.expected_prezero_ms=1000*(T.stim_trigger_getsecs-T.stim_software_start_getsecs);
T.audio_t1_actual_minus_scheduled_ms=1000*(T.actual_pattern_onset-T.scheduled_pattern_onset);
T.audio_t5_actual_minus_scheduled_ms=1000*(T.actual_fifth_tone_onset-T.scheduled_fifth_tone_onset);
T.t5_ttl_write_vs_actual_t5_ms=1000*(T.t5_ttl_write_getsecs-T.actual_fifth_tone_onset);
T.expected_stim_vs_t5_ttl_write_ms=1000*(T.stim_trigger_getsecs-T.t5_ttl_write_getsecs);
T.saved_offset_residual_ms=T.stim_t5_offset_ms-T.expected_stim_vs_scheduled_t5_ms;
T.saved_actual_offset_residual_ms=T.stim_actual_t5_offset_ms-T.expected_stim_vs_actual_t5_ms;
stim=T.is_stimulated==1;
T.stim_timing_valid=stim & isfinite(T.expected_stim_vs_scheduled_t5_ms);
T.stim_outside_tolerance=T.stim_timing_valid & ...
    abs(T.expected_stim_vs_scheduled_t5_ms-cfg.stimTargetMs)>cfg.toleranceMs;
T.stim_near_minus400=T.stim_timing_valid & abs(T.expected_stim_vs_scheduled_t5_ms+400)<=50;
S=table();
allMetrics={'audio_t1_actual_minus_scheduled_ms','audio_t5_actual_minus_scheduled_ms','t5_ttl_write_vs_actual_t5_ms'};
for j=1:numel(allMetrics)
    S=[S;erta_timing_stats(T.(allMetrics{j}),allMetrics{j},'all log trials',0,cfg.toleranceMs)]; %#ok<AGROW>
end
metrics={'expected_stim_vs_scheduled_t5_ms','expected_stim_vs_actual_t5_ms', ...
    'software_start_vs_t5_ms','expected_prezero_ms','expected_stim_vs_t5_ttl_write_ms', ...
    'saved_offset_residual_ms','saved_actual_offset_residual_ms'};
targets=[cfg.stimTargetMs cfg.stimTargetMs NaN NaN NaN 0 0];
if isfield(info,'stim') && isfield(info.stim,'preZeroSec')
    targets(3)=cfg.stimTargetMs-1000*info.stim.preZeroSec;
    targets(4)=1000*info.stim.preZeroSec;
end
for j=1:numel(metrics)
    S=[S;erta_timing_stats(T.(metrics{j})(stim),metrics{j},'stimulated log trials',targets(j),cfg.toleranceMs)]; %#ok<AGROW>
end
for b=unique(T.block(stim & isfinite(T.block)))'
    S=[S;erta_timing_stats(T.expected_stim_vs_scheduled_t5_ms(stim & T.block==b), ...
        'expected_stim_vs_scheduled_t5_ms',sprintf('stimulated block %g',b),cfg.stimTargetMs,cfg.toleranceMs)]; %#ok<AGROW>
end
F.stim_absent=zeros(height(F),1); F.stim_NaN=F.stim_absent;
F.stim_Inf=F.stim_absent; F.stim_invalid_or_empty=F.stim_absent;
for j=1:height(F)
    [v,a,b]=erta_numeric_field(E,char(F.field(j)));
    F.stim_absent(j)=sum(a&stim); F.stim_NaN(j)=sum(isnan(v)&~a&~b&stim);
    F.stim_Inf(j)=sum(isinf(v)&stim); F.stim_invalid_or_empty(j)=sum(b&stim);
end
fprintf('\nACQUISITION LOG: %d trials; %d stimulated; %d unknown stimulation flags.\n',N,nnz(stim),nnz(~ismember(T.is_stimulated,[0 1])));
fprintf('Expected onset near -400 ms (+/-50): %d; valid stimulated timings: %d/%d.\n', ...
    nnz(T.stim_near_minus400),nnz(T.stim_timing_valid),nnz(stim));
fprintf('Offsets beyond target +/- %.3f ms: %d. Positive offsets mean LATE.\n',cfg.toleranceMs,nnz(T.stim_outside_tolerance));
valid=T.expected_stim_vs_scheduled_t5_ms(T.stim_timing_valid);
if ~isempty(valid)
    fprintf('Expected stimulation vs scheduled T5: mean %+.3f ms, median %+.3f ms, SD %.3f ms, range [%+.3f, %+.3f] ms.\n', ...
        mean(valid),median(valid),std(valid,0),min(valid),max(valid));
end
fprintf('These are log-derived estimates, not measurements of delivered current.\n');
nonSoftware=stim & (ismissing(T.stim_delivery_mode) | ~strcmpi(T.stim_delivery_mode,'software'));
if any(nonSoftware) || any(stim & T.stim_hardware_available~=1)
    warning('Some stimulated rows are not confirmed software-mode/hardware-available rows. Inspect per-trial mode and hardware flags; do not interpret simulation timing as delivery.');
end
disp(S);
end

function [v,absent,invalid]=erta_numeric_field(E,name)
v=nan(numel(E),1); absent=repmat(~isfield(E,name),numel(E),1); invalid=false(size(v));
if all(absent), return; end
for j=1:numel(E)
    z=E(j).(name);
    if (isnumeric(z)||islogical(z)) && isscalar(z) && isreal(z), v(j)=double(z);
    else, invalid(j)=true; end
end
end

function S=erta_timing_stats(x,metric,scope,target,tolerance)
x=double(x(:)); f=sort(x(isfinite(x))); n=numel(f); q=nan(1,7);
if n>0
    % Linear quantiles, no Statistics Toolbox needed.
    p=[0 .05 .25 .5 .75 .95 1]; ix=1+(n-1)*p;
    q=f(floor(ix))'.*(ceil(ix)-ix)+f(ceil(ix))'.*(ix-floor(ix));
    exact=floor(ix)==ceil(ix); q(exact)=f(ix(exact));
end
avg=NaN; sd=NaN; outside=NaN;
if n>0, avg=mean(f); if n>1, sd=std(f,0); end; end
if isfinite(target), outside=sum(abs(f-target)>tolerance); end
S=table(string(metric),string(scope),numel(x),n,sum(isnan(x)),sum(isinf(x)), ...
    avg,sd,q(1),q(2),q(3),q(4),q(5),q(6),q(7),target,tolerance,outside, ...
    'VariableNames',{'metric','scope','n_total','n_valid','n_nan','n_inf','mean_ms','sd_ms', ...
    'min_ms','p05_ms','p25_ms','median_ms','p75_ms','p95_ms','max_ms','target_ms','tolerance_ms','n_outside_tolerance'});
end

function [Q,S]=erta_tone_marker_timing(Q,code,pos,fs,cfg)
S=table(); if isempty(Q), return; end
for tone=1:4
    x=nan(height(Q),1); counts=zeros(height(Q),1);
    for j=1:height(Q)
        lower=-inf; if j>1, lower=Q.t5_sample(j-1); end
        candidates=pos(code==tone+4 & pos>lower & pos<Q.t5_sample(j));
        counts(j)=numel(candidates);
        if numel(candidates)==1, x(j)=candidates; end
    end
    Q.(sprintf('t%d_marker_count',tone))=counts;
    Q.(sprintf('t%d_sample',tone))=x;
    offset=1000*(x-Q.t5_sample)/fs;
    Q.(sprintf('t%d_relative_t5_ms',tone))=offset;
    err=offset-cfg.toneOnsetsMs(tone);
    Q.(sprintf('t%d_timing_error_ms',tone))=err;
    S=[S;erta_timing_stats(err,sprintf('T%d relative to T5: error from expected interval',tone), ...
        'recorded EEG markers',0,cfg.toleranceMs)]; %#ok<AGROW>
end
S=[S;erta_timing_stats(Q.condition_offset_ms,'S9 to condition LABEL','recorded EEG markers',NaN,cfg.toleranceMs); ...
    erta_timing_stats(Q.stim_offset_ms(Q.stim_paired),'S9 to S10 LABEL (not current onset)', ...
    'recorded EEG markers',NaN,cfg.toleranceMs)];
Q.s10_label_offset_ms=Q.stim_offset_ms;
end
