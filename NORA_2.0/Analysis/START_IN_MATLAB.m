function plan = START_IN_MATLAB(overrides)
% Portable, sequential ERtaVNS master. Developed for MATLAB R2025a.
% Put this file in NORA_2.0 or NORA_2.0/Analysis, then press Run.
% All participants finish a stage before the next stage begins.
% Analysis scripts and raw recordings are not modified. Existing analysis
% outputs may be overwritten by the analysis scripts, as in a normal rerun.

%% ================= EDIT THESE SETTINGS ===============================
cfg.scriptsToRun = { ...
    'A1_ERtaVNS_marker_QC'
    'B2_EEG_preparation_LG'
    'B3_EEG_preprocessing_LG'
    'C1_Pupil_preprocessing'
    'E0_ERtaVNS_trial_count_waveform_stability'
}; % Reorder/remove lines; add other self-contained, noninteractive scripts.

cfg.projectRoot = '';  % Auto: this file's project, then MATLAB Current Folder.
% Optional explicit root: 'L:\Downloads\ERtaVNS\NORA_2.0'
%                     or 'F:\ERtaVNS\NORA_2.0'
cfg.eeglabRoot = '';   % Auto-detect; or enter the folder containing eeglab.m.
cfg.showFigures = false; % Figures are still saved by the analysis scripts.
cfg.dryRun = false;      % true: display/check the queue without running it.
%% ====================================================================

% Optional function-call overrides are useful for testing without editing.
if nargin>0
    assert(isstruct(overrides) && isscalar(overrides),'Overrides must be a scalar struct.');
    names=fieldnames(overrides);
    for k=1:numel(names)
        assert(isfield(cfg,names{k}),'Unknown master setting: %s',names{k});
        cfg.(names{k})=overrides.(names{k});
    end
end
root=find_project(cfg.projectRoot,fileparts(mfilename('fullpath')),pwd);
analysis=fullfile(root,'Analysis'); raw=fullfile(root,'Experiments','Data');
scripts=cellstr(string(cfg.scriptsToRun)); scripts=scripts(:);
assert(~isempty(scripts),'scriptsToRun is empty.');
for i=1:numel(scripts)
    scripts{i}=regexprep(scripts{i},'\.m$','');
    assert(~isempty(regexp(scripts{i},'^[A-Za-z]\w*$','once')), ...
        'Use a script name, without a folder: %s',scripts{i});
    assert(isfile(fullfile(analysis,[scripts{i} '.m'])),'Missing script: %s',scripts{i});
end
[records,sessions]=recordings(raw);
visits=unique([sessions.visit]);
needsEEG=any(ismember(scripts,{'B2_EEG_preparation_LG','B3_EEG_preprocessing_LG'}));
eegRoot=''; lookup='';
if needsEEG, [eegRoot,lookup]=find_eeglab(cfg.eeglabRoot,root); end
plan=make_plan(scripts,records,sessions);
assert(~isempty(plan),'No runnable jobs were found.');
% Prepare and validate every known script before starting any analysis.
sources=cell(numel(plan),1);
for i=1:numel(plan)
    original=fileread(fullfile(analysis,[plan(i).script '.m']));
    sources{i}=batch_source(original,plan(i),root,analysis,eegRoot,lookup,visits);
    check_inputs(plan(i),records,sessions,original);
end
fprintf('\nPROJECT: %s\nOUTPUTS: %s\n',root,fullfile(root,'Analysis_Output'));
if needsEEG, fprintf('EEGLAB: %s\nElectrode locations: %s\n',eegRoot,lookup); end
for i=1:numel(plan), fprintf('%2d. %s\n',i,plan(i).label); end
if cfg.dryRun
    fprintf('\nDRY RUN: queue and batch settings checked; no analyses or files written.\n');
    return;
end

out=fullfile(root,'Analysis_Output');
runDir=fullfile(out,'Master_Run_Logs',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
assert(~isfolder(runDir),'This run folder already exists: %s',runDir);
[ok,msg]=mkdir(runDir); assert(ok,'Cannot create run log folder: %s',msg);
oldPath=path; oldFolder=pwd; oldVisible=get(groot,'DefaultFigureVisible');
restore=onCleanup(@()restore_environment(oldPath,oldFolder,oldVisible)); %#ok<NASGU>
if needsEEG, select_eeglab(eegRoot); end
addpath(analysis,'-begin');
if cfg.showFigures, set(groot,'DefaultFigureVisible','on');
else, set(groot,'DefaultFigureVisible','off'); end
state=struct('status','RUNNING','project',root,'job',0,'total',numel(plan), ...
    'current','','updated','','error','');
save_state(runDir,state,plan);
for i=1:numel(plan)
    state.job=i; state.current=plan(i).label;
    plan(i).status='RUNNING'; plan(i).started=datestr(now,31);
    save_state(runDir,state,plan);
    diary(fullfile(runDir,sprintf('%02d_%s.log',i,plan(i).script)));
    fprintf('\nSTART %d/%d: %s | %s\n',i,numel(plan),state.current,plan(i).started);
    try
        if known_script(plan(i).script)
            jobDir=fullfile(runDir,sprintf('batch_%02d',i)); mkdir(jobDir);
            scriptFile=fullfile(jobDir,[plan(i).script '.m']);
            write_text(scriptFile,sources{i});
            % An unexpected missing file must raise an error, not open a
            % selection dialog and silently wait overnight.
            write_text(fullfile(jobDir,'erta_no_dialog.m'),sprintf([ ...
                'function varargout=erta_no_dialog(varargin)\n' ...
                'error(''ERtaVNS:BatchInput'',''Missing batch input. Inspect the run log; file selection is disabled.'');\nend\n']));
        else
            scriptFile=fullfile(analysis,[plan(i).script '.m']);
        end
        started=now;
        execute_isolated(scriptFile);
        verify_products(plan(i),records,sessions,out,started,sources{i});
        plan(i).status='COMPLETE'; plan(i).finished=datestr(now,31);
        fprintf('\nCOMPLETE: %s\n',state.current);
        diary off;
        save_state(runDir,state,plan);
    catch ME
        plan(i).status='FAILED'; plan(i).finished=datestr(now,31);
        state.status='FAILED'; state.error=getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'\nSTOPPED: %s\n%s\n',state.current,state.error);
        save_state(runDir,state,plan); diary off;
        fprintf(2,'Run log: %s\n',runDir);
        rethrow(ME);
    end
end
state.status='COMPLETE'; state.current='All selected stages';
save_state(runDir,state,plan);
fprintf('\nALL SELECTED STAGES COMPLETE.\nOutputs: %s\nLogs: %s\n',out,runDir);
end

function execute_isolated(scriptFile)
% Synchronous execution: run returns only after the script finishes.
% Its clear/clc commands cannot erase the master function's variables.
run(scriptFile);
end

function root=find_project(explicit,here,current)
if ~isempty(explicit)
    root=char(explicit); assert(is_project(root),'Invalid projectRoot: %s',root); return;
end
for start={here,current}
    p=start{1};
    while ~isempty(p)
        if is_project(p), root=p; return; end
        parent=fileparts(p); if strcmp(parent,p), break; end; p=parent;
    end
end
error('ERtaVNS:Project','Put this master in NORA_2.0 or Analysis, or set cfg.projectRoot.');
end
function yes=is_project(p)
yes=isfolder(fullfile(p,'Analysis')) && isfolder(fullfile(p,'Experiments','Data'));
end

function [R,S]=recordings(raw)
R=struct('participant',{},'visit',{},'task',{},'folder',{},'stem',{});
S=struct('participant',{},'visit',{},'folder',{});
subs=dir(raw); subs=subs([subs.isdir]); [~,ix]=sort({subs.name}); subs=subs(ix);
for si=1:numel(subs)
    p=subs(si).name; if isempty(regexp(p,'^S\d+$','once')), continue; end
    vv=dir(fullfile(raw,p,'Visit*')); vv=vv([vv.isdir]);
    for vi=1:numel(vv)
        tok=regexp(vv(vi).name,'^Visit([12])$','tokens','once'); if isempty(tok), continue; end
        v=str2double(tok{1}); folder=fullfile(raw,p,vv(vi).name);
        S(end+1)=struct('participant',p,'visit',v,'folder',folder); %#ok<AGROW>
        for task={'E1','E2'}
            stem=sprintf('%sV%d%s',p,v,task{1});
            if any(cellfun(@(ext)isfile(fullfile(folder,[stem ext])),{'.vmrk','.vhdr','.asc'}))
                R(end+1)=struct('participant',p,'visit',v,'task',task{1},'folder',folder,'stem',stem); %#ok<AGROW>
            end
        end
    end
end
end

function jobs=make_plan(scripts,R,S)
blank=struct('script','','label','','visit',0,'participant','','task','', ...
    'status','PENDING','started','','finished',''); jobs=repmat(blank,0,1);
for k=1:numel(scripts)
    name=scripts{k}; j=blank; j.script=name; j.label=name;
    switch name
        case 'A1_ERtaVNS_marker_QC'
            assert(~isempty(R),'No E1/E2 recordings were found for A1.');
            for r=1:numel(R)
                j.visit=R(r).visit; j.participant=R(r).participant; j.task=R(r).task;
                j.label=[name ' | ' R(r).stem]; jobs(end+1)=j; %#ok<AGROW>
            end
        case {'B2_EEG_preparation_LG','B3_EEG_preprocessing_LG'}
            assert(~isempty(R),'No E1/E2 recordings were found for %s.',name);
            for v=unique([S.visit])
                if ~any([R.visit]==v), continue; end
                j.visit=v; j.label=sprintf('%s | Visit%d | all participants',name,v);
                jobs(end+1)=j; %#ok<AGROW>
            end
        otherwise
            jobs(end+1)=j; %#ok<AGROW>
    end
end
end

function yes=known_script(name)
yes=ismember(name,{'A1_ERtaVNS_marker_QC','B2_EEG_preparation_LG', ...
    'B3_EEG_preprocessing_LG','C1_Pupil_preprocessing','E0_ERtaVNS_trial_count_waveform_stability'});
end

function s=batch_source(s,j,root,analysis,eegRoot,lookup,visits)
if ~known_script(j.script), return; end
if startsWith(j.script,'E0_')
    s=assign(s,'cfg.projectRoot',quoted(root));
else
    if startsWith(j.script,'A1_'), s=assign(s,'THIS_DIR',quoted(analysis));
    else, s=assign(s,'LG_THIS_SCRIPT_DIR',quoted(analysis)); end
    s=assign(s,'PROJECT_ROOT',quoted(root));
end
switch j.script
    case 'A1_ERtaVNS_marker_QC'
        replacement=sprintf('answ = {%s,%s,%s};',quoted(j.participant),quoted(num2str(j.visit)),quoted(j.task));
        s=replace_once(s,'(?s)answ\s*=\s*inputdlg\(.*?\);',replacement,'A1 recording dialog');
        s=assign(s,'acquisitionMatOverride','''''');
        s=strrep(s,'uigetfile(','erta_no_dialog(');
    case {'B2_EEG_preparation_LG','B3_EEG_preprocessing_LG'}
        s=assign(s,'analysisVisit',num2str(j.visit));
        s=assign(s,'mainpath',quoted([eegRoot filesep]));
        if startsWith(j.script,'B2_'), s=assign(s,'startSubIndex','1');
        else, s=assign(s,'num_subject','1'); end
        pat='\[mainpath\s*,?\s*''plugins[\\/]dipfit[^'']*standard-10-5-cap385\.elp''\]';
        [a,b]=regexp(s,pat,'start','end');
        assert(~isempty(a),'Unsupported electrode-location expression in %s.',j.script);
        for k=numel(a):-1:1, s=[s(1:a(k)-1) quoted(lookup) s(b(k)+1:end)]; end
    case 'C1_Pupil_preprocessing'
        assert(~isempty(visits),'No Visit1/Visit2 folders were found for C1.');
        s=assign(s,'visitsToAnalyze',mat2str(visits));
        s=replace_once(s,'(?m)^startSubIdx\s*=\s*\d+\s*;', ...
            'startSubIdx = 1;','C1 initial participant index');
        s=replace_once(s,'(?m)^endSubIdx\s*=\s*(inf|\d+)\s*;', ...
            'endSubIdx = inf;','C1 final participant index');
        pat='disp\(getReport\(ME,''extended''\)\);\s*continue;';
        s=replace_once(s,pat,'rethrow(ME);','C1 task error handler');
end
% Preserve scientific/QC warnings, but make caught processing errors fatal.
ids={'B2:ControlFailed','B3:ControlFailed','B3:TaskMetadata','B3:NewFigure', ...
    'B3:ICLabelFailed','C1:ControlSkipped','E0:Skipped','E0:FigureFailed'};
for k=1:numel(ids), s=strrep(s,['warning(''' ids{k} ''''],['error(''' ids{k} '''']); end
end
function s=assign(s,name,value)
pat=['(?m)^' regexptranslate('escape',name) '\s*=\s*[^;\r\n]*;'];
s=replace_once(s,pat,[name ' = ' value ';'],[name ' setting']);
end
function s=replace_once(s,pat,replacement,description)
[a,b]=regexp(s,pat,'start','end');
assert(numel(a)==1,'ERtaVNS:ScriptVersion', ...
    'Cannot safely adapt %s (%d matches). This analysis version needs a batch adapter update.',description,numel(a));
s=[s(1:a-1) replacement s(b+1:end)];
end
function q=quoted(s)
q=['''' strrep(char(s),'''','''''') ''''];
end

function [root,lookup]=find_eeglab(explicit,project)
candidates={};
if ~isempty(explicit), candidates={char(explicit)};
else
    w=which('eeglab'); if ~isempty(w), candidates{end+1}=fileparts(w); end
    user=getenv('USERPROFILE'); if isempty(user), user=getenv('HOME'); end
    parents={fullfile(project,'Programs'),fileparts(project),fullfile(user,'Desktop'),fullfile(user,'Downloads')};
    drive=regexp(project,'^[A-Za-z]:','match','once');
    if ~isempty(drive), parents{end+1}=fullfile([drive filesep],'Programs'); end
    for p=1:numel(parents)
        d=dir(fullfile(parents{p},'eeglab*')); d=d([d.isdir]);
        for k=1:numel(d), candidates{end+1}=fullfile(d(k).folder,d(k).name); end %#ok<AGROW>
    end
end
root=''; lookup='';
for k=1:numel(candidates)
    r=candidates{k}; if ~isfile(fullfile(r,'eeglab.m')), continue; end
    preferred=fullfile(r,'plugins','dipfit5.4','standard_BESA','standard-10-5-cap385.elp');
    if isfile(preferred), lookup=preferred;
    else
        d=dir(fullfile(r,'plugins','dipfit*','standard_BESA','standard-10-5-cap385.elp'));
        if ~isempty(d), lookup=fullfile(d(1).folder,d(1).name); end
    end
    if ~isempty(lookup), root=r; return; end
end
error('ERtaVNS:EEGLAB','Set cfg.eeglabRoot to an EEGLAB folder with the DIPFIT electrode-location file.');
end
function select_eeglab(root)
other=which('eeglab','-all'); if ischar(other), other={other}; end
entries=strsplit(path,pathsep);
for k=1:numel(other)
    p=fileparts(other{k}); if strcmpi(p,root), continue; end
    for e=1:numel(entries)
        if strcmpi(entries{e},p) || startsWith(lower(entries{e}),lower([p filesep]))
            rmpath(entries{e});
        end
    end
end
addpath(root,'-begin');
end

function check_inputs(j,R,S,source)
if startsWith(j.script,'A1_')
    rows=find(strcmp({R.participant},j.participant) & [R.visit]==j.visit & strcmp({R.task},j.task));
    for r=rows, require_files(R(r).folder,R(r).stem,{'.vmrk','.vhdr','.mat'}); end
elseif startsWith(j.script,'B2_') || startsWith(j.script,'B3_')
    for r=find([R.visit]==j.visit)
        require_files(R(r).folder,R(r).stem,{'.vhdr','.vmrk','.eeg'});
        if startsWith(j.script,'B3_'), require_files(R(r).folder,R(r).stem,{'.mat'}); end
    end
    controlFlag='prepareControl'; if startsWith(j.script,'B3_'), controlFlag='processControl'; end
    if setting_enabled(source,controlFlag)
        for si=find([S.visit]==j.visit)
            stem=sprintf('%sV%dC',S(si).participant,S(si).visit);
            if isfile(fullfile(S(si).folder,[stem '.vhdr']))
                require_files(S(si).folder,stem,{'.vhdr','.vmrk','.eeg'});
            end
        end
    end
elseif strcmp(j.script,'C1_Pupil_preprocessing')
    assert(~isempty(R),'No E1/E2 recordings were found for C1.');
    for r=1:numel(R), require_files(R(r).folder,R(r).stem,{'.asc'}); end
end
end
function yes=setting_enabled(source,name)
tok=regexp(source,['(?m)^' name '\s*=\s*(true|false|1|0)\s*;'],'tokens','once');
assert(~isempty(tok),'Cannot determine batch setting %s; use true/false in the analysis script.',name);
yes=ismember(tok{1},{'true','1'});
end
function require_files(folder,stem,exts)
for k=1:numel(exts)
    f=fullfile(folder,[stem exts{k}]); d=dir(f);
    assert(isscalar(d) && d.bytes>0,'Missing or empty input: %s',f);
end
end

function verify_products(j,R,S,out,t0,source)
% Check fresh core products; caught plot/ICLabel failures are fatal above.
switch j.script
    case 'A1_ERtaVNS_marker_QC'
        stem=sprintf('%sV%d%s',j.participant,j.visit,j.task);
        fresh(fullfile(out,'Pilot_QC',[stem '_marker_QC.csv']),t0);
        if contains(source,'_onset_timing_report.txt')
            fresh(fullfile(out,'Pilot_QC',[stem '_onset_timing_report.txt']),t0);
        end
    case 'B2_EEG_preparation_LG'
        if setting_enabled(source,'prepareTasks')
            rr=R([R.visit]==j.visit);
            for p=unique({rr.participant})
                fresh(fullfile(out,sprintf('Visit%d',j.visit),p{1},['ica_' p{1} '.set']),t0);
            end
        end
        if setting_enabled(source,'prepareControl'), verify_controls(j,S,out,t0,'ica_','_C.set','.vhdr'); end
    case 'B3_EEG_preprocessing_LG'
        for r=find([R.visit]==j.visit)
            folder=fullfile(out,sprintf('Visit%d',j.visit),R(r).participant);
            fresh(fullfile(folder,sprintf('auditory_all_LG_%s_%s.set',R(r).participant,R(r).task)),t0);
        end
        if setting_enabled(source,'processControl'), verify_controls(j,S,out,t0,'control_all_','.set','.vhdr'); end
    case 'C1_Pupil_preprocessing'
        for r=1:numel(R)
            fresh(fullfile(out,sprintf('Visit%d',R(r).visit),R(r).participant, ...
                sprintf('%s_V%d_%s_pupil_preproc.mat',R(r).participant,R(r).visit,R(r).task)),t0);
        end
        for si=1:numel(S)
            stem=sprintf('%sV%dC',S(si).participant,S(si).visit);
            if isfile(fullfile(S(si).folder,[stem '.asc']))
                fresh(fullfile(out,sprintf('Visit%d',S(si).visit),S(si).participant, ...
                    sprintf('%s_V%d_C_pupil_preproc.mat',S(si).participant,S(si).visit)),t0);
            end
        end
    case 'E0_ERtaVNS_trial_count_waveform_stability'
        logs=dir(fullfile(out,'Visit*','S*','Trial_Count_Stability','run_log.csv'));
        assert(~isempty(logs),'No stability run logs were created.');
        for k=1:numel(logs)
            f=fullfile(logs(k).folder,logs(k).name);
            if logs(k).datenum<t0-2/86400, continue; end % unrelated old dataset
            T=readtable(f,'TextType','string');
            assert(~isempty(T) && all(T.status=="OK"),'Stability run contains missing/failed modalities: %s',f);
        end
        assert(any([logs.datenum]>=t0-2/86400),'No fresh stability logs were created.');
end
end
function verify_controls(j,S,out,t0,prefix,suffix,rawExt)
for si=find([S.visit]==j.visit)
    stem=sprintf('%sV%dC',S(si).participant,S(si).visit);
    if isfile(fullfile(S(si).folder,[stem rawExt]))
        fresh(fullfile(out,sprintf('Visit%d',j.visit),S(si).participant,[prefix S(si).participant suffix]),t0);
    end
end
end
function fresh(file,t0)
d=dir(file); assert(isscalar(d) && d.bytes>0,'Missing/empty output: %s',file);
assert(d.datenum>=t0-2/86400,'Output was not refreshed in this job: %s',file);
end
function save_state(folder,state,plan)
state.updated=datestr(now,31);
write_text(fullfile(folder,'run_status.json'),jsonencode(struct('run',state,'jobs',plan),'PrettyPrint',true));
write_text(fullfile(folder,'RUN_STATUS.txt'),sprintf('%s\nJob %d/%d: %s\nUpdated: %s\n%s\n', ...
    state.status,state.job,state.total,state.current,state.updated,state.error));
save(fullfile(folder,'run_plan.mat'),'plan','state');
end
function write_text(file,content)
f=fopen(file,'w','n','UTF-8'); assert(f>=0,'Cannot write: %s',file);
c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s',content);
end
function restore_environment(oldPath,oldFolder,oldVisible)
diary off; path(oldPath); cd(oldFolder); set(groot,'DefaultFigureVisible',oldVisible);
end
