%% E0_ERtaVNS_trial_count_waveform_stability
% Put in NORA_2.0/Analysis and run after B3/C1. MATLAB R2020b+.
% Includes ALL retained trials, including habituation; excludes buffer flags.
% No re-filtering, re-baselining, S10 alignment or changes to raw data.
% Scans every available participant, Visit1/Visit2 and E1/E2 automatically.
%
% Per dataset: GD and GS waveform figures, GD and GS projection figures,
% plus separate GD/GS correlation and RMSE figures. PNG ONLY.
% All files go in Analysis_Output/Visit#/Sxx/Trial_Count_Stability/.
% No group directory, no FIG files, no extra task subfolders are created.
%
% n = clean trials averaged. FULL = every retained non-buffer trial.
% Shading = middle 95% of repeated random subset averages for the smallest n
% displayed in that panel (pointwise resampling interval, NOT SEM).
% Metric plots summarize 1000 random samples WITHOUT replacement per n.
% Waveform display curves use one reproducible nested draw, as in old E0.
%
% Shorter-experiment projection: keep ALL available habituation fixed and
% draw floor(clean TEST trials * N/180) test trials, N=100:10:180 per block.
% This assumes unchanged test-trial yield and exchangeable test trials.
% Exact planned critical count ACROSS TWO BLOCKS = 0.2*N per GD category
% (also paired GS). These planned counts precede artifact rejection.
% For the mandatory-buffer design, planned pure GS = N test + 50 habituation.
% Legacy projections follow the observed legacy mix; they do not validate
% the current mandatory-buffer design. Never derive N+1 from cleaned rows.
% A subset overlaps FULL: agreement is a stability check, not independent
% reliability/power. No automatic acceptance cutoff is imposed.
%
% Optional override: E0_config=struct('projectRoot','X:\...\NORA_2.0');
% The script leaves existing output files in place and never deletes them.

cfg=struct;
cfg.projectRoot=fileparts(fileparts(mfilename('fullpath')));
cfg.nResamples=1000;
cfg.seed=20260925;
cfg.retainedCounts=[10 15 20 25 30];
cfg.candidateTrialsPerBlock=100:10:180;
cfg.referenceTrialsPerBlock=180;
cfg.erpChannels={'Pz'};
cfg.pupilPlotWindow=[-500 2000];
cfg.erpPlotWindow=[-1400 700]; % Whole EEG epoch: -1.4 to +0.7 seconds (axis in ms).
cfg.pupilMetricWindow=[0 2000]; cfg.erpMetricWindow=[0 700];
cfg.saveFigures=true; cfg.figureVisible='off';
% ================= PLOT SCALE: CHANGE HERE =============================
cfg.pupilUpperLimit=1.2; % Change this number to your preferred pupil maximum.
cfg.pupilYLim=[];       % [] = automatic shared lower edge + maximum above.
                       % Set e.g. [-0.5 1.2] to fix BOTH pupil limits instead.
cfg.erpYLim=[];         % [] = automatic shared ERP scale; or set [min max].
% A numeric pupilYLim overrides pupilUpperLimit. These display settings do
% not change the data or the post-T5 windows used for correlation/RMSE.
% ======================================================================
if exist('E0_config','var')
    assert(isstruct(E0_config) && isscalar(E0_config),'E0_config must be a scalar struct.');
    names=fieldnames(E0_config);
    for k=1:numel(names)
        assert(isfield(cfg,names{k}),'Unknown E0_config setting: %s',names{k});
        cfg.(names{k})=E0_config.(names{k});
    end
end
assert(cfg.nResamples>=2 && cfg.nResamples==fix(cfg.nResamples),'nResamples must be an integer >=2.');
assert(cfg.referenceTrialsPerBlock==180,'Reference design is 180 test trials per block.');
assert(all(cfg.candidateTrialsPerBlock>=100 & cfg.candidateTrialsPerBlock<=180 & mod(cfg.candidateTrialsPerBlock,10)==0),'Use block lengths 100:10:180.');
assert(all(isfinite(cfg.retainedCounts) & cfg.retainedCounts>=1 & cfg.retainedCounts==fix(cfg.retainedCounts)),'Invalid clean trial counts.');
assert(isscalar(cfg.seed) && isfinite(cfg.seed) && cfg.seed>=0 && cfg.seed==fix(cfg.seed),'Invalid random seed.');
for name={'pupilPlotWindow','erpPlotWindow','pupilMetricWindow','erpMetricWindow'}
    w=cfg.(name{1}); assert(numel(w)==2 && all(isfinite(w)) && w(2)>w(1),'Invalid time window.');
end
assert(cfg.pupilMetricWindow(1)>=0 && cfg.erpMetricWindow(1)>=0,'Metrics use post-T5 windows.');
for name={'pupilYLim','erpYLim'}
    w=cfg.(name{1}); assert(isempty(w) || (numel(w)==2 && all(isfinite(w)) && w(2)>w(1)),'Invalid y-axis limits.');
end
assert(isscalar(cfg.pupilUpperLimit) && isfinite(cfg.pupilUpperLimit),'Invalid pupil upper limit.');
analysisRoot=fullfile(cfg.projectRoot,'Analysis_Output');
datasets=discover_datasets(analysisRoot);
assert(~isempty(datasets),'No preprocessed E1/E2 products found. Put E0 in Analysis or set E0_config.projectRoot.');
runLog=table; allSummary=table; allAudit=table;
savedFiles=strings(0,1); observedLimits=[Inf -Inf; Inf -Inf];
for di=1:numel(datasets)
    ds=datasets(di); outDir=fullfile(ds.folder,'Trial_Count_Stability');
    if ~isfolder(outDir), mkdir(outDir); end
    fprintf('\n%s | %s | Visit %d\n',ds.participant,ds.task,ds.visit);
    try
        [condition,conditionSource]=resolve_condition(ds,cfg);
        cells=cell(2,1); audits=cell(2,1);
        for mi=1:2
            mods=["Pupil","ERP"]; modality=mods(mi);
            try
                if mi==1, [cells{mi},audits{mi}]=load_pupil(ds,condition,cfg);
                else, [cells{mi},audits{mi}]=load_erp(ds,cfg); end
                runLog=[runLog;log_row(ds,modality,"OK",audits{mi}.buffer_status)]; %#ok<AGROW>
            catch ME
                runLog=[runLog;log_row(ds,modality,"SKIPPED",string(ME.message))]; %#ok<AGROW>
                warning('E0:Skipped','%s %s: %s',ds.participant,modality,ME.message);
            end
        end
        if ~isempty(audits{1}) && ~isempty(audits{2}) && audits{2}.buffer_status=="legacy_or_rule_violation"
            audits{1}.buffer_status="legacy_or_rule_violation";
        end
        results=cell(2,2,2); datasetSummary=table; datasetAudit=table;
        for mi=1:2
            if isempty(cells{mi}), continue; end
            C=cells{mi}; A=audits{mi};
            A.participant=ds.participant; A.visit=ds.visit; A.task=ds.task; A.session_condition=condition;
            datasetAudit=[datasetAudit;struct2table(A)]; %#ok<AGROW>
            for gi=1:2
                for si=1:2
                    gs=["GD","GS"]; sc=["pure","paired"];
                    mask=C.global==gs(gi) & C.stim==(si==2) & C.keep;
                    key=ds.participant+"|"+ds.visit+"|"+ds.task+"|"+condition+"|"+C.modality+"|"+gs(gi)+"|"+sc(si);
                    R=analyse_cell(C.X(mask,:),C.hab(mask),C.t,C.metricWindow,C.plotWindow,key,gs(gi),sc(si),cfg);
                    R.modality=C.modality; R.units=C.units; R.buffer_status=A.buffer_status;
                    results{mi,gi,si}=R; T=R.summary;
                    T.participant=repmat(ds.participant,height(T),1); T.visit=repmat(ds.visit,height(T),1);
                    T.task=repmat(ds.task,height(T),1); T.session_condition=repmat(condition,height(T),1);
                    T.condition_source=repmat(conditionSource,height(T),1); T.modality=repmat(C.modality,height(T),1);
                    T.global_status=repmat(gs(gi),height(T),1); T.stim_class=repmat(sc(si),height(T),1);
                    T.buffer_status=repmat(A.buffer_status,height(T),1); T.source_file=repmat(C.source,height(T),1);
                    T.projection_basis=repmat("observed_clean_test_yield_with_fixed_habituation",height(T),1);
                    if A.buffer_status=="legacy_or_rule_violation"
                        T.projection_basis(:)="observed_legacy_trial_mix_with_fixed_habituation";
                    end
                    datasetSummary=[datasetSummary;T]; %#ok<AGROW>
                end
            end
        end
        if ~isempty(datasetSummary)
            stem=sprintf('%s_V%d_%s_%s',ds.participant,ds.visit,ds.task,condition);
            matFile=fullfile(outDir,[stem '_results.mat']);
            writetable(datasetSummary,fullfile(outDir,[stem '_metrics.csv']));
            save(matFile,'results','datasetSummary','datasetAudit','cfg','ds','condition','conditionSource','-v7.3');
            savedFiles(end+1)=string(matFile); %#ok<SAGROW>
            lim=waveform_limits(results);
            observedLimits(:,1)=min(observedLimits(:,1),lim(:,1));
            observedLimits(:,2)=max(observedLimits(:,2),lim(:,2));
            allSummary=[allSummary;datasetSummary]; allAudit=[allAudit;datasetAudit]; %#ok<AGROW>
        end
    catch ME
        runLog=[runLog;log_row(ds,"Dataset","SKIPPED",string(ME.message))]; %#ok<AGROW>
        warning('E0:Skipped','%s %s: %s',ds.participant,ds.task,ME.message);
    end
    writetable(runLog(runLog.participant==ds.participant & runLog.visit==ds.visit,:),fullfile(outDir,'run_log.csv'));
end
% A second pass gives every participant/task figure exactly the same pupil
% scale and the same ERP scale, without retaining all waveforms in memory.
plotLimits=observedLimits;
for mi=1:2
    if any(~isfinite(plotLimits(mi,:))), plotLimits(mi,:)=[-1 1]; end
    span=diff(plotLimits(mi,:)); if span==0, span=1; end
    plotLimits(mi,:)=plotLimits(mi,:)+[-1 1]*.06*span;
end
if isempty(cfg.pupilYLim)
    plotLimits(1,2)=cfg.pupilUpperLimit;
    if plotLimits(1,1)>=plotLimits(1,2), plotLimits(1,1)=plotLimits(1,2)-1; end
else
    plotLimits(1,:)=cfg.pupilYLim;
end
if ~isempty(cfg.erpYLim), plotLimits(2,:)=cfg.erpYLim; end
if cfg.saveFigures
    for fi=1:numel(savedFiles)
        S=load(savedFiles(fi),'results','ds','condition');
        try
            dataset_figures(S.results,S.ds,S.condition,fileparts(savedFiles(fi)),cfg,plotLimits);
            save(savedFiles(fi),'plotLimits','-append');
        catch ME
            runLog=[runLog;log_row(S.ds,"Figures","FAILED",string(ME.message))]; %#ok<AGROW>
            writetable(runLog(runLog.participant==S.ds.participant & runLog.visit==S.ds.visit,:),fullfile(fileparts(savedFiles(fi)),'run_log.csv'));
            warning('E0:FigureFailed','%s',ME.message);
        end
    end
end
fprintf('\nFinished: %d datasets. PNGs and metrics are in each participant''s Trial_Count_Stability folder.\n',numel(savedFiles));

%% Local functions
function D=discover_datasets(root)
D=struct('participant',{},'visit',{},'task',{},'folder',{},'pupilFile',{},'eegFile',{});
vv=dir(fullfile(root,'Visit*'));
for v=1:numel(vv)
    tok=regexp(vv(v).name,'^Visit([12])$','tokens','once');
    if ~vv(v).isdir || isempty(tok), continue; end
    visit=str2double(tok{1}); pp=dir(fullfile(root,vv(v).name,'S*'));
    for p=1:numel(pp)
        if ~pp(p).isdir || isempty(regexp(pp(p).name,'^S\d+$','once')), continue; end
        folder=fullfile(root,vv(v).name,pp(p).name);
        for task=["E1","E2"]
            pf=fullfile(folder,sprintf('%s_V%d_%s_pupil_preproc.mat',pp(p).name,visit,task));
            ef=fullfile(folder,sprintf('auditory_all_LG_%s_%s.set',pp(p).name,task));
            if ~isfile(pf) && ~isfile(ef), continue; end
            D(end+1)=struct('participant',string(pp(p).name),'visit',visit,'task',task, ...
                'folder',folder,'pupilFile',pf,'eegFile',ef); %#ok<AGROW>
        end
    end
end
end

function [condition,source]=resolve_condition(ds,cfg)
values=strings(0,1); sources=strings(0,1);
if isfile(ds.pupilFile)
    P=load(ds.pupilFile,'dataTableClean','analysisVisit');
    if isfield(P,'analysisVisit'), assert(isequal(double(P.analysisVisit),ds.visit),'Pupil analysisVisit/path mismatch.'); end
    if isfield(P,'dataTableClean') && istable(P.dataTableClean)
        T=P.dataTableClean;
        if ismember('visit',T.Properties.VariableNames)
            x=double(T.visit); assert(all(~isfinite(x) | x==ds.visit),'Pupil visit/path mismatch.');
        end
        if ismember('session_condition',T.Properties.VariableNames)
            x=normal_condition(T.session_condition);
            values=[values; x]; sources=[sources; repmat("pupil_metadata",numel(x),1)];
        end
    end
end
% Read actual saved tracker rows; never invoke assignment generation.
fp=fullfile(cfg.projectRoot,'Experiments','ERtaVNS_participant_tracker.csv');
if isfile(fp)
    T=readtable(fp,'TextType','string'); field=sprintf('Visit%dCondition',ds.visit);
    if all(ismember({'Participant',field},T.Properties.VariableNames))
        rows=strcmpi(strtrim(string(T.Participant)),ds.participant);
        assert(nnz(rows)<=1,'Duplicate participant tracker rows.');
        x=normal_condition(T.(field)(rows)); values=[values;x];
        sources=[sources;repmat("tracker",numel(x),1)];
    end
end
u=unique(values);
assert(numel(u)==1,'ACTIVE/SHAM metadata missing or conflicting; dataset skipped rather than guessed.');
condition=u(1); source=strjoin(unique(sources,'stable'),"+");
end

function x=normal_condition(v)
x=upper(strtrim(string(v(:))));
x=x(~ismissing(x) & x~="" & x~="NAN" & x~="<UNDEFINED>");
assert(all(ismember(x,["ACTIVE","SHAM"])),'Unrecognized session_condition metadata.');
end

function [C,A]=load_pupil(ds,condition,cfg)
assert(isfile(ds.pupilFile),'Pupil file missing.');
P=load(ds.pupilFile,'dataTableClean','relTimes','removedTrialsSummaryAll');
assert(isfield(P,'dataTableClean') && isfield(P,'relTimes'),'Missing C1 output variables.');
T=P.dataTableClean; need={'effect','stim_paired','PUPIL_Z','modality','is_post_stim_buffer'};
assert(istable(T) && all(ismember(need,T.Properties.VariableNames)),'Missing pupil labels, phase or buffer flags; rerun current C1.');
t=double(P.relTimes(:)'); X=cells_matrix(T.PUPIL_Z,numel(t));
stim=binary_values(T.stim_paired); buf=binary_values(T.is_post_stim_buffer);
phase=lower(strtrim(string(T.modality))); effect=effect_labels(T.effect);
assert(all(ismember(phase,["test","habituation"])),'Unknown pupil phase.');
if ismember('session_condition',T.Properties.VariableNames)
    x=normal_condition(T.session_condition); assert(all(x==condition),'Pupil condition conflict.');
end
keep=~buf;
A=audit_template("Pupil",height(T),nnz(phase=="habituation"),nnz(buf));
A.buffer_status="flags_present_rule_unverified";
% C1 saves exclusions but does not save its complete raw buffer-rule report.
% Exclusion evidence proves buffers were handled, not that the entire N+1
% rule was verified. Record this provenance with the per-dataset metrics.
if isfield(P,'removedTrialsSummaryAll') && istable(P.removedTrialsSummaryAll) && ...
        ismember('Reason',P.removedTrialsSummaryAll.Properties.VariableNames)
    A.n_buffers_removed_upstream=nnz(contains(lower(string(P.removedTrialsSummaryAll.Reason)),"post-stimulation"));
    if A.n_buffers_removed_upstream>0, A.buffer_status="buffer_exclusions_recorded"; end
end
if any(buf), A.buffer_status="buffer_exclusions_recorded"; end
if ismember('trialInFile',T.Properties.VariableNames)
    assert(numel(unique(T.trialInFile))==height(T),'Duplicate pupil trial identities.');
end
C=make_cell(X,t,effect,stim,keep,phase=="habituation","Pupil","z",string(ds.pupilFile),cfg.pupilMetricWindow,cfg.pupilPlotWindow);
A.n_eligible=nnz(keep); A.n_habituation_included=nnz(C.hab & keep);
end

function [C,A]=load_erp(ds,cfg)
assert(isfile(ds.eegFile),'EEG .set missing.');
% Read standard EEGLAB MAT .set directly. Supports embedded numeric data and
% external little-endian float32 .fdt; avoids changing EEGLAB paths/state.
S=load(ds.eegFile,'-mat');
if isfield(S,'EEG'), E=S.EEG; else, E=S; end
assert(isfield(E,'epoch') && numel(E.epoch)==E.trials && E.trials>0,'No EEG epoch metadata.');
need={'eventlg_effect','eventlg_stim_paired','eventlg_phase','eventlg_post_stim_buffer','eventtype','eventlatency'};
assert(all(isfield(E.epoch,need)),'Missing scalar B3 epoch metadata; rerun current B3.');
effect=strings(E.trials,1); stim=false(E.trials,1); buf=stim; phase=effect;
for k=1:E.trials
    ep=E.epoch(k);
    effect(k)=scalar_text(ep.eventlg_effect); phase(k)=upper(scalar_text(ep.eventlg_phase));
    stim(k)=scalar_binary(ep.eventlg_stim_paired); buf(k)=scalar_binary(ep.eventlg_post_stim_buffer);
    if isfield(ep,'eventlg_exclude_primary'), buf(k)=buf(k) | scalar_binary(ep.eventlg_exclude_primary); end
    types=upper(strtrim(string(ep.eventtype))); lats=numeric_list(ep.eventlatency);
    assert(numel(types)==numel(lats),'Inconsistent epoch event types/latencies.');
    near=ismember(types(:),["FLS","FLD","RLS","RLD"]) & abs(lats(:))<=max(1,1.5*1000/E.srate);
    assert(nnz(near)==1,'Missing/ambiguous zero-latency T5 LG anchor; rerun B3.');
    if isfield(ep,'eventoriginal_type')
        original=string(ep.eventoriginal_type); original=original(:);
        assert(numel(original)==numel(near) && ~isempty(regexp(char(original(near)),'^S\s*9$','once')), ...
            'The zero-latency LG event is not an original S9 marker.');
    end
end
assert(all(ismember(phase,["TEST","HABITUATION"])),'Missing/unknown EEG phase; no habituation guessing.');
labels=string({E.chanlocs.labels}); ci=zeros(1,numel(cfg.erpChannels));
for k=1:numel(ci)
    ix=find(strcmpi(labels,cfg.erpChannels{k}));
    assert(numel(ix)==1,'Missing/duplicate ERP channel: %s',cfg.erpChannels{k}); ci(k)=ix;
end
if ~isnumeric(E.data)
    dataName=char(string(E.data)); [~,name,ext]=fileparts(dataName);
    assert(strcmpi(ext,'.fdt'),'Unsupported external EEG format; expected .fdt.');
    fp=fullfile(ds.folder,[name ext]);
    fid=fopen(fp,'r','ieee-le'); assert(fid>=0,'Missing EEG .fdt: %s',fp);
    closer=onCleanup(@()fclose(fid));
    [data,nread]=fread(fid,Inf,'single=>single'); clear closer;
    assert(nread==double(E.nbchan)*double(E.pnts)*double(E.trials),'Unexpected .fdt length.');
    E.data=reshape(data,[E.nbchan E.pnts E.trials]);
end
assert(numel(E.data)==double(E.nbchan)*double(E.pnts)*double(E.trials),'EEG data dimensions disagree.');
X=reshape(mean(double(E.data(ci,:,:)),1,'omitnan'),E.pnts,E.trials).';
keep=~buf;
A=audit_template("ERP",E.trials,nnz(phase=="HABITUATION"),nnz(buf));
A.buffer_status="flags_present_rule_unverified";
if isfield(E,'etc') && isfield(E.etc,'LG_trigger_report')
    Q=E.etc.LG_trigger_report;
    if isfield(Q,'bufferRuleDetected')
        if scalar_binary(Q.bufferRuleDetected), A.buffer_status="complete_rule_detected_upstream";
        else, A.buffer_status="legacy_or_rule_violation"; end
    end
    if isfield(Q,'nPostStimBuffers'), A.n_buffers_detected_upstream=double(Q.nPostStimBuffers); end
end
if A.buffer_status=="flags_present_rule_unverified" && any(buf), A.buffer_status="buffer_exclusions_recorded"; end
A.n_eligible=nnz(keep); A.n_habituation_included=nnz(phase=="HABITUATION" & keep);
C=make_cell(X,double(E.times(:)'),effect_labels(effect),stim,keep,phase=="HABITUATION","ERP","uV",string(ds.eegFile),cfg.erpMetricWindow,cfg.erpPlotWindow);
end

function A=audit_template(modality,n,hab,buf)
A=struct('modality',modality,'n_input',n,'n_habituation_input',hab,'n_habituation_included',0,'n_habituation_excluded',0,'n_buffer_flags_excluded',buf, ...
    'n_buffers_removed_upstream',NaN,'n_buffers_detected_upstream',NaN,'n_eligible',0,'buffer_status',"");
end
function C=make_cell(X,t,effect,stim,keep,hab,modality,units,source,metricWindow,plotWindow)
assert(size(X,2)==numel(t) && size(X,1)==numel(keep),'Waveform/metadata dimensions disagree.');
assert(all(isfinite(t)) && all(diff(t)>0),'Invalid time axis.');
% EEGLAB half-open epochs can end one sample before the requested endpoint.
dt=median(diff(t));
assert(t(1)<=metricWindow(1)+dt && t(end)>=metricWindow(2)-dt,'Data do not cover the requested metric window (within one sample).');
assert(nnz(t>=metricWindow(1) & t<=metricWindow(2))>=3,'Fewer than three metric-window samples.');
assert(~any(isinf(X(:))),'Infinite waveform sample; inspect preprocessing.');
C=struct('X',X,'t',t,'global',effect(:),'stim',stim(:),'keep',keep(:),'hab',hab(:), ...
    'modality',modality,'units',units,'source',source,'metricWindow',metricWindow,'plotWindow',plotWindow);
end
function s=scalar_text(v)
while iscell(v) && numel(v)==1, v=v{1}; end
s=strtrim(string(v)); assert(isscalar(s) && ~ismissing(s),'Expected scalar metadata.');
end
function b=binary_values(v)
if iscell(v), v=cellfun(@(x)double(binary_values(x)),v); end
if islogical(v), b=v(:); return; end
if ~isnumeric(v)
    v=lower(strtrim(string(v))); ok=ismember(v,["0","1","false","true"]);
    assert(all(ok(:)),'Missing or invalid boolean metadata.'); b=v=="1" | v=="true";
else
    assert(all(isfinite(v(:)) & (v(:)==0 | v(:)==1)),'Missing or invalid boolean metadata.'); b=logical(v);
end
b=b(:);
end
function b=scalar_binary(v)
b=binary_values(v); assert(isscalar(b),'Missing/ambiguous scalar boolean metadata.');
end
function x=numeric_list(v)
if iscell(v), x=cellfun(@(z)str2double(string(z)),v); else, x=double(v); end
end
function g=effect_labels(v)
v=lower(strtrim(string(v))); g=strings(size(v));
g(v=="global standard")="GS"; g(v=="global deviant")="GD";
assert(all(g~=""),'Unknown global-status labels.');
end
function X=cells_matrix(v,nt)
if isempty(v), X=zeros(0,nt); return; end
if iscell(v)
    assert(all(cellfun(@numel,v)==nt),'Pupil sample counts differ.');
    X=zeros(numel(v),nt); for k=1:numel(v), X(k,:)=double(v{k}(:)'); end
else, X=double(v); end
end

function R=analyse_cell(X,hab,t,win,plotwin,key,gs,sc,cfg)
n=size(X,1); nHab=nnz(hab); nTest=nnz(~hab); metricMask=t>=win(1) & t<=win(2);
R=struct('time_ms',t,'full',mean(X,1,'omitnan'),'n_available',n,'n_habituation',nHab, ...
    'metric_window',win,'plot_window',plotwin,'draws',[],'summary',table,'key',key);
% A clean-count row draws from all retained trials. A projection row keeps
% habituation fixed, reducing ONLY the test phase of the same recording.
P=table;
for count=cfg.retainedCounts
    P=[P;table("clean_count",NaN,NaN,count,-1,NaN,'VariableNames', ...
        {'analysis','trials_per_block','planned_critical_n','n','n_habituation','n_test'})]; %#ok<AGROW>
end
for N=cfg.candidateTrialsPerBlock
    cleanTest=floor(nTest*N/cfg.referenceTrialsPerBlock);
    P=[P;table("projection",N,.2*N,nHab+cleanTest,nHab,cleanTest,'VariableNames',P.Properties.VariableNames)]; %#ok<AGROW>
end
% Key by both total count and fixed habituation count: equal total counts
% need not have the same composition. -1 means a simple random clean subset.
spec=unique([P.n P.n_habituation],'rows');
D=struct('n',{},'n_habituation',{},'seed',{},'r',{},'rmse',{},'example',{},'lo_curve',{},'hi_curve',{});
displayStream=RandStream('mt19937ar','Seed',key_seed(key+"|display",cfg.seed));
order=randperm(displayStream,n); hids=find(hab); tids=find(~hab); testOrder=tids(randperm(displayStream,nTest));
for k=1:size(spec,1)
    count=spec(k,1); fixedHab=spec(k,2);
    if count<1 || count>n, continue; end
    seed=key_seed(key+"|n="+count+"|hab="+fixedHab,cfg.seed);
    stream=RandStream('mt19937ar','Seed',seed); ids=zeros(cfg.nResamples,count);
    for b=1:cfg.nResamples
        if fixedHab<0, ids(b,:)=randperm(stream,n,count);
        else, ids(b,:)=[hids(:)' tids(randperm(stream,nTest,count-fixedHab))']; end
    end
    W=sparse(repmat((1:cfg.nResamples)',count,1),ids(:),1,cfg.nResamples,n);
    finite=isfinite(X); X0=X; X0(~finite)=0;
    denom=W*double(finite); Y=(W*X0)./denom; Y(denom==0)=NaN;
    ref=R.full(metricMask); [r,rmse]=wave_metrics(Y(:,metricMask),ref,isfinite(ref));
    bands=column_quantiles(Y,[.025 .975]);
    if fixedHab<0, chosen=order(1:count);
    else, chosen=[hids(:)' testOrder(1:count-fixedHab)']; end
    D(end+1)=struct('n',count,'n_habituation',fixedHab,'seed',seed,'r',r,'rmse',rmse, ...
        'example',mean(X(chosen,:),1,'omitnan'),'lo_curve',bands(1,:),'hi_curve',bands(2,:)); %#ok<AGROW>
end
for k=1:height(P)
    qr=nan(1,3); qe=qr; done=0; nr=0; ne=0; seed=NaN; drawIndex=NaN; status="ok";
    if P.n(k)>n, status="insufficient_clean_trials";
    elseif P.n(k)<1, status="no_clean_trials";
    elseif nnz(isfinite(R.full(metricMask)))<3, status="insufficient_samples";
    else
        drawIndex=find([D.n]==P.n(k) & [D.n_habituation]==P.n_habituation(k));
        d=D(drawIndex); done=cfg.nResamples; nr=nnz(isfinite(d.r)); ne=nnz(isfinite(d.rmse));
        qr=finite_quantile(d.r,[.025 .5 .975]); qe=finite_quantile(d.rmse,[.025 .5 .975]); seed=d.seed;
        if nr<done || ne<done, status="some_metrics_undefined"; end
    end
    row=P(k,:); row.n_available=n; row.n_available_test=nTest; row.n_available_habituation=nHab;
    row.n_resamples_completed=done; row.n_valid_r=nr; row.n_valid_RMSE=ne;
    row.median_r=qr(2); row.r_low95=qr(1); row.r_high95=qr(3);
    row.median_RMSE=qe(2); row.RMSE_low95=qe(1); row.RMSE_high95=qe(3);
    row.window_start_ms=win(1); row.window_end_ms=win(2);
    samples=t(metricMask); row.actual_window_start_ms=min(samples); row.actual_window_end_ms=max(samples);
    row.seed=seed; row.status=status; row.draw_index=drawIndex;
    row.planned_test_cell_n=NaN; row.planned_habituation_cell_n=NaN;
    if row.analysis=="projection"
        if gs=="GS" && sc=="pure"
            row.planned_test_cell_n=row.trials_per_block; row.planned_habituation_cell_n=50;
        else
            row.planned_test_cell_n=row.planned_critical_n; row.planned_habituation_cell_n=0;
        end
    end
    R.summary=[R.summary;row];
end
R.draws=D;
end
function [r,e,nvalid]=wave_metrics(Y,ref,support)
% Fixed reference support across repeats; no changing time-window pairs.
Y=Y(:,support); ref=ref(support); nvalid=sum(isfinite(Y),2);
r=nan(size(Y,1),1); e=r;
if numel(ref)<3, return; end
ok=all(isfinite(Y),2); Z=Y(ok,:);
e(ok)=sqrt(mean((Z-ref).^2,2));
a=ref-mean(ref); z=Z-mean(Z,2); den=sqrt(sum(a.^2)*sum(z.^2,2));
v=sum(z.*a,2)./den; v(den==0)=NaN;
finite=isfinite(v); v(finite)=max(-1,min(1,v(finite))); r(ok)=v;
end
function q=column_quantiles(X,p)
q=nan(numel(p),size(X,2));
if all(isfinite(X(:)))
    X=sort(X,1); h=1+(size(X,1)-1)*p(:); lo=floor(h); hi=ceil(h);
    q=X(lo,:).*(1-(h-lo))+X(hi,:).*(h-lo);
else
    for j=1:size(X,2), q(:,j)=finite_quantile(X(:,j),p).'; end
end
end
function q=finite_quantile(x,p)
x=sort(x(isfinite(x))); q=nan(size(p)); if isempty(x), return; end
h=1+(numel(x)-1)*p; lo=floor(h); hi=ceil(h);
q=reshape(x(lo),size(p)).*(1-(h-lo))+reshape(x(hi),size(p)).*(h-lo);
end
function seed=key_seed(key,base)
% Stable hash in exactly representable integer range; independent of loop order.
h=mod(double(base),2147483646);
for b=double(char(key)), h=mod(h*131+b,2147483646); end
seed=h+1;
end

function lim=waveform_limits(results)
lim=[Inf -Inf; Inf -Inf];
for mi=1:2
    for gi=1:2
        for si=1:2
            R=results{mi,gi,si}; if isempty(R), continue; end
            mask=R.time_ms>=R.plot_window(1) & R.time_ms<=R.plot_window(2);
            values=R.full(mask);
            for k=1:numel(R.draws), values=[values R.draws(k).example(mask)]; end %#ok<AGROW>
            for mode=["clean_count","projection"]
                T=R.summary(R.summary.analysis==mode & R.summary.n_resamples_completed>0,:);
                if isempty(T), continue; end
                [~,j]=min(T.n); d=R.draws(T.draw_index(j));
                values=[values d.lo_curve(mask) d.hi_curve(mask)]; %#ok<AGROW>
            end
            values=values(isfinite(values));
            if ~isempty(values), lim(mi,:)=[min(lim(mi,1),min(values)) max(lim(mi,2),max(values))]; end
        end
    end
end
end
function dataset_figures(results,ds,condition,out,cfg,plotLimits)
stem=sprintf('%s_V%d_%s_%s',ds.participant,ds.visit,ds.task,condition);
for gi=1:2
    globals=["GD","GS"]; g=globals(gi);
    for mode=["clean_count","projection"]
        f=figure('Color','w','Visible',cfg.figureVisible,'Position',[60 50 1450 950]); cleaner=onCleanup(@()close(f));
        tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
        for mi=1:2
            for si=1:2
                ax=nexttile(tl,(mi-1)*2+si); R=results{mi,gi,si}; ttl=panel_title(mi,si,g,condition,cfg);
                if isempty(R)
                    title(ax,ttl); text(ax,.1,.5,'Data unavailable','Units','normalized'); ylim(ax,plotLimits(mi,:));
                    if mi==1, xlim(ax,cfg.pupilPlotWindow); else, xlim(ax,cfg.erpPlotWindow); end
                    continue;
                end
                waveform_panel(ax,R,mode,ttl,plotLimits(mi,:));
            end
        end
        label='Trial-count waveform stability';
        if mode=="projection", label='Trial-count waveform stability | 100-180 test trials/block'; end
        title(tl,sprintf('%s | %s | %s | Visit %d | %s | %s | T5 = 0 ms',ds.participant,ds.task,condition,ds.visit,g,label),'Interpreter','none');
        suffix='waveforms'; if mode=="projection", suffix='projections'; end
        exportgraphics(f,fullfile(out,sprintf('%s_%s_%s.png',stem,g,suffix)),'Resolution',180); clear cleaner;
    end
    for metric=1:2
        f=figure('Color','w','Visible',cfg.figureVisible,'Position',[60 50 1450 950]); cleaner=onCleanup(@()close(f));
        tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
        % Same metric scale for pure/paired within each modality.
        metricLims=[-1 1.02;-1 1.02];
        if metric==2
            for mi=1:2
                upper=0;
                for si=1:2
                    R=results{mi,gi,si}; if isempty(R), continue; end
                    v=R.summary.RMSE_high95(R.summary.analysis=="projection"); v=v(isfinite(v));
                    if ~isempty(v), upper=max(upper,max(v)); end
                end
                if upper==0, upper=1; end
                metricLims(mi,:)=[0 upper*1.1];
            end
        end
        for mi=1:2
            for si=1:2
                ax=nexttile(tl,(mi-1)*2+si); R=results{mi,gi,si}; title(ax,panel_title(mi,si,g,condition,cfg));
                if isempty(R), text(ax,.1,.5,'Data unavailable','Units','normalized'); continue; end
                T=sortrows(R.summary(R.summary.analysis=="projection",:),'trials_per_block');
                if metric==1, y=T.median_r; lo=T.r_low95; hi=T.r_high95; ylabel(ax,'Correlation with FULL (r)');
                else, y=T.median_RMSE; lo=T.RMSE_low95; hi=T.RMSE_high95; ylabel(ax,"RMSE ("+R.units+")"); end
                errorbar(ax,T.trials_per_block,y,y-lo,hi-y,'o-','LineWidth',1.5);
                xlabel(ax,'Test trials per block'); grid(ax,'on'); xticks(ax,cfg.candidateTrialsPerBlock);
                xlim(ax,[min(cfg.candidateTrialsPerBlock)-3 max(cfg.candidateTrialsPerBlock)+3]); ylim(ax,metricLims(mi,:));
                text(ax,.03,.08,'Dots: median | bars: middle 95% of repeated samples','Units','normalized', ...
                    'FontSize',8,'BackgroundColor','w','VerticalAlignment','bottom','Margin',2,'Clipping','on');
            end
        end
        labels=["Correlation with FULL","RMSE compared with FULL"]; files=["correlation","RMSE"];
        title(tl,sprintf('%s | %s | %s | Visit %d | %s | %s | %d random samples',ds.participant,ds.task,condition,ds.visit,g,labels(metric),cfg.nResamples),'Interpreter','none');
        exportgraphics(f,fullfile(out,sprintf('%s_%s_%s.png',stem,g,files(metric))),'Resolution',180); clear cleaner;
    end
end
end
function ttl=panel_title(mi,si,g,condition,cfg)
mods=["Pupil","ERP"]; classes=["Pure",condition+"-paired"];
ttl=mods(mi)+" | "+classes(si)+" "+g;
if mi==2, ttl=ttl+" | "+strjoin(string(cfg.erpChannels),', '); end
end
function waveform_panel(ax,R,mode,ttl,limits)
hold(ax,'on'); T=R.summary(R.summary.analysis==mode,:);
T=T(T.n_resamples_completed>0 & T.n<R.n_available,:);
nCurves=max(1,height(T)); colors=lines(nCurves); note="";
% Avoid repeating MATLAB's seven default colors when more curves are shown.
if nCurves>7
    colors=hsv2rgb([(0:nCurves-1)'/nCurves repmat(.75,nCurves,1) repmat(.8,nCurves,1)]);
end
if ~isempty(T)
    [~,j]=min(T.n); d=R.draws(T.draw_index(j));
    fill(ax,[R.time_ms fliplr(R.time_ms)],[d.lo_curve fliplr(d.hi_curve)],[.72 .8 .9], ...
        'EdgeColor','none','FaceAlpha',.35,'HandleVisibility','off');
    note=sprintf('Shading: middle 95%% of random n=%d averages',d.n);
    for j=1:height(T)
        d=R.draws(T.draw_index(j)); label=sprintf('n=%d',d.n);
        if mode=="projection", label=sprintf('%d trials/block | n=%d',T.trials_per_block(j),d.n); end
        plot(ax,R.time_ms,d.example,'Color',colors(j,:),'LineWidth',1.25,'DisplayName',label);
    end
end
fullLabel=sprintf('FULL n=%d',R.n_available);
if mode=="projection", fullLabel=sprintf('180 trials/block | FULL n=%d',R.n_available); end
plot(ax,R.time_ms,R.full,'Color',[0 0 0],'LineWidth',3,'DisplayName',fullLabel);
xline(ax,0,'--','T5','HandleVisibility','off'); yline(ax,0,':','HandleVisibility','off');
xlim(ax,R.plot_window); ylim(ax,limits); grid(ax,'on'); xlabel(ax,'Time from T5 onset (ms)');
set(ax,'Clipping','on');
if R.modality=="Pupil", yticks(ax,unique([yticks(ax) limits(2)])); end
if R.modality=="Pupil", ylabel(ax,'Pupil (z)'); else, ylabel(ax,'\muV'); end
title(ax,sprintf('%s | FULL n=%d',ttl,R.n_available));
legend(ax,'Location','best','Box','off','FontSize',8,'Interpreter','none');
if R.n_available==0, note="No retained trials"; end
if strlength(note)>0
    text(ax,.02,.08,note,'Units','normalized','FontSize',8,'BackgroundColor','w','Interpreter','none', ...
        'VerticalAlignment','bottom','Margin',2,'Clipping','on');
end
end
function T=log_row(ds,modality,status,message)
T=table(ds.participant,ds.visit,ds.task,modality,status,message, ...
    'VariableNames',{'participant','visit','task','modality','status','message'});
end
