function out = parse_LG_eyelink(lines, varargin)
% PARSE_LG_EYELINK
% Parse CURRENT ER-taVNS EyeLink messages from A1_Exp_01_VNS / B1_Exp_02_VNS.
%
% Tone 1-4 example:
%   MSG <time> LG_TONE TRIAL_001 BLOCK_1 POS_1 CODE_5 PHYS_A
%
% T5 examples:
%   E1: MSG <time> LG_T5 ... CODE_9 PHASE_TEST TYPE_FLS LOCAL_LS GLOBAL_GS
%       PHYS_A PATTERN_AAAAA STIM_1 SESSION_ACTIVE VISIT_1
%   E2: MSG <time> LG_T5 ... CODE_9 PHASE_TEST TYPE_FLS LOCAL_LS GLOBAL_GS
%       PHYS_A PATTERN_AAAAA SESSION_ACTIVE VISIT_1
%
% On stimulation-paired trials both tasks also send:
%   MSG <time> VNS_TRIAL TRIAL_001 BLOCK_1 CODE_10 CONDITION_ACTIVE SITE_...
%
% The parser therefore derives stimulation status from either the optional
% STIM_ field on LG_T5 or the VNS_TRIAL message.  No 50-ms timing correction
% is applied: messages are written at the intended tone onset.

p=inputParser;
p.addParameter('expectedTrials',410,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('verbose',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
expectedTrials=p.Results.expectedTrials;
verbose=p.Results.verbose;

if ischar(lines)||isstring(lines), lines=cellstr(lines); end

pat14=['^MSG\s+([0-9]+(?:\.[0-9]+)?)\s+LG_TONE\s+' ...
       'TRIAL_(\d+)\s+BLOCK_(\d+)\s+POS_([1-4])\s+' ...
       'CODE_(\d+)\s+PHYS_([AB])\s*$'];

% Optional STIM_n followed by mandatory SESSION/VISIT.
pat5=['^MSG\s+([0-9]+(?:\.[0-9]+)?)\s+LG_T5\s+' ...
      'TRIAL_(\d+)\s+BLOCK_(\d+)\s+POS_5\s+' ...
      'CODE_(\d+)\s+PHASE_([A-Z]+)\s+TYPE_([A-Z_]+)\s+' ...
      'LOCAL_([A-Z]+)\s+GLOBAL_([A-Z]+)\s+' ...
      'PHYS_([AB])\s+PATTERN_([AB]{5})' ...
      '(?:\s+STIM_(\d+))?\s+SESSION_([A-Z]+)\s+VISIT_(\d+)\s*$'];


patBuffer=['^MSG\s+([0-9]+(?:\.[0-9]+)?)\s+POST_STIM_BUFFER\s+' ...
        'TRIAL_(\d+)\s+BLOCK_(\d+)\s+CODE_(\d+)\s+TYPE_([A-Z]+)\s+' ...
        'EXCLUDE_EEG_1\s+EXCLUDE_PUPIL_1\s*$'];

patVNS=['^MSG\s+([0-9]+(?:\.[0-9]+)?)\s+VNS_TRIAL\s+' ...
        'TRIAL_(\d+)\s+BLOCK_(\d+)\s+CODE_(\d+)\s+' ...
        'CONDITION_([A-Z]+)\s+SITE_(.+)\s*$'];

parsed=struct('time',{},'trial',{},'block',{},'pos',{},'code',{},'phys',{}, ...
    'phase',{},'rawType',{},'local',{},'global',{},'pattern',{}, ...
    'stim',{},'session',{},'visit',{});
bufferTrials=[]; bufferTimes=[];
vnsTrials=[]; vnsTimes=[]; vnsConditions=strings(0,1); vnsSites=strings(0,1);
maxTrial=0;

for i=1:numel(lines)
    s=strtrim(lines{i});
    tok=regexp(s,pat14,'tokens','once');
    if ~isempty(tok)
        e.time=str2double(tok{1}); e.trial=str2double(tok{2}); e.block=str2double(tok{3});
        e.pos=str2double(tok{4}); e.code=str2double(tok{5}); e.phys=tok{6};
        e.phase=''; e.rawType=''; e.local=''; e.global=''; e.pattern='';
        e.stim=NaN; e.session=''; e.visit=NaN;
        parsed(end+1)=e; %#ok<AGROW>
        maxTrial=max(maxTrial,e.trial); continue;
    end

    tok=regexp(s,pat5,'tokens','once');
    if ~isempty(tok)
        e.time=str2double(tok{1}); e.trial=str2double(tok{2}); e.block=str2double(tok{3});
        e.pos=5; e.code=str2double(tok{4}); e.phase=tok{5}; e.rawType=tok{6};
        e.local=tok{7}; e.global=tok{8}; e.phys=tok{9}; e.pattern=tok{10};
        % MATLAB regexp omits the unmatched optional STIM_(n) token entirely.
        % Therefore LG_T5 messages without STIM_ contain 12 tokens, whereas
        % messages with STIM_ contain 13 tokens.
        if numel(tok) >= 13
            e.stim=str2double(tok{11});
            e.session=tok{12};
            e.visit=str2double(tok{13});
        else
            e.stim=NaN;
            e.session=tok{11};
            e.visit=str2double(tok{12});
        end
        parsed(end+1)=e; %#ok<AGROW>
        maxTrial=max(maxTrial,e.trial); continue;
    end

    tok=regexp(s,patBuffer,'tokens','once');
    if ~isempty(tok)
        bufferTimes(end+1,1)=str2double(tok{1}); %#ok<AGROW>
        bufferTrials(end+1,1)=str2double(tok{2}); %#ok<AGROW>
        maxTrial=max(maxTrial,str2double(tok{2})); continue;
    end

    tok=regexp(s,patVNS,'tokens','once');
    if ~isempty(tok)
        vnsTimes(end+1,1)=str2double(tok{1}); %#ok<AGROW>
        vnsTrials(end+1,1)=str2double(tok{2}); %#ok<AGROW>
        vnsConditions(end+1,1)=string(tok{5}); %#ok<AGROW>
        vnsSites(end+1,1)=string(tok{6}); %#ok<AGROW>
        maxTrial=max(maxTrial,str2double(tok{2}));
    end
end

nTrials=max(maxTrial,expectedTrials);
tmsg=nan(nTrials,5); code=nan(nTrials,5); phys=strings(nTrials,5);
block=nan(nTrials,1); phase=strings(nTrials,1); rawType=strings(nTrials,1);
local=strings(nTrials,1); globalLab=strings(nTrials,1); pattern=strings(nTrials,1);
t5ms=nan(nTrials,1); t5code=nan(nTrials,1); stimFlag=false(nTrials,1);
sessionCondition=strings(nTrials,1); visit=nan(nTrials,1); stimSite=strings(nTrials,1);
postStimBuffer=false(nTrials,1);
duplicateSlots=0;

for i=1:numel(parsed)
    e=parsed(i);
    if e.trial<1||e.trial>nTrials||e.pos<1||e.pos>5, continue; end
    if isfinite(tmsg(e.trial,e.pos)), duplicateSlots=duplicateSlots+1; end
    tmsg(e.trial,e.pos)=e.time; code(e.trial,e.pos)=e.code; phys(e.trial,e.pos)=string(e.phys); block(e.trial)=e.block;
    if e.pos==5
        phase(e.trial)=string(e.phase); rawType(e.trial)=string(e.rawType);
        local(e.trial)=string(e.local); globalLab(e.trial)=string(e.global);
        pattern(e.trial)=string(e.pattern); t5ms(e.trial)=e.time; t5code(e.trial)=e.code;
        if isfinite(e.stim), stimFlag(e.trial)=logical(e.stim); end
        sessionCondition(e.trial)=string(e.session); visit(e.trial)=e.visit;
    end
end

% VNS_TRIAL is authoritative for paired-stimulation status in both tasks.
for i=1:numel(vnsTrials)
    tr=vnsTrials(i);
    if tr>=1&&tr<=nTrials
        stimFlag(tr)=true;
        if strlength(sessionCondition(tr))==0, sessionCondition(tr)=vnsConditions(i); end
        stimSite(tr)=vnsSites(i);
    end
end

% Buffer identity is derived from the experimental sequence, NOT from a
% dedicated POST_STIM_BUFFER message.  For every stimulation-paired trial N,
% N+1 must exist, be a GS, and be unstimulated.  We apply N+1 exclusion only
% when the complete recording obeys this rule.  This protects legacy pilots
% (e.g. recordings acquired before the mandatory buffer rule).
stimIdx = find(stimFlag);
bufferCandidate = false(nTrials,1);
bufferRuleViolation = false(nTrials,1);

for k=1:numel(stimIdx)
    tr = stimIdx(k);
    nx = tr + 1;
    if nx > nTrials
        bufferRuleViolation(tr)=true;
        continue;
    end

    sameBlock = isfinite(block(tr)) && isfinite(block(nx)) && block(nx)==block(tr);
    nextIsGS = upper(globalLab(nx))=="GS";
    nextUnstim = ~stimFlag(nx);

    if sameBlock && nextIsGS && nextUnstim
        bufferCandidate(nx)=true;
    else
        bufferRuleViolation(tr)=true;
    end
end

bufferRuleDetected = ~isempty(stimIdx) && ~any(bufferRuleViolation(stimIdx));
if bufferRuleDetected
    postStimBuffer = bufferCandidate;
else
    postStimBuffer(:)=false;
end

% Historical POST_STIM_BUFFER messages, if present, are retained only for QC.
explicitBufferMessage = false(nTrials,1);
for i=1:numel(bufferTrials)
    tr=bufferTrials(i);
    if tr>=1 && tr<=nTrials, explicitBufferMessage(tr)=true; end
end

% Current acquisition already writes FLS/FLD/RLS/RLD at T5 for both TEST
% and HABITUATION.  Phase remains available separately.
typeShort=upper(rawType);
trial=(1:nTrials).';
isMatchedGS = stimFlag & upper(globalLab)=="GS";
isStimGD    = stimFlag & upper(globalLab)=="GD";

trialInfo=table(trial,block,phase,rawType,typeShort,local,globalLab,pattern, ...
    phys(:,5),t5code,t5ms,stimFlag,isMatchedGS,isStimGD,postStimBuffer,postStimBuffer,postStimBuffer,sessionCondition,visit,stimSite, ...
    'VariableNames',{'trial','block','phase','raw_type_short','type_short','local_short', ...
    'global_short','physical_pattern','fifth_tone_physical','t5_code','t5_time_ms', ...
    'stim_paired','is_matched_gs_control','is_stimulated_gd','is_post_stim_buffer','exclude_primary_eeg','exclude_primary_pupil','session_condition','visit','stim_site'});

report=struct();
report.expectedTrials=expectedTrials; report.nTrialsIndexed=maxTrial;
report.nParsedToneMessages=numel(parsed); report.nVNSMessages=numel(vnsTrials);
report.duplicateToneSlots=duplicateSlots; report.completeToneTrials=sum(all(isfinite(tmsg),2));
report.missingToneTrials=find(~all(isfinite(tmsg),2)); report.nT5=sum(isfinite(t5ms));
report.nTest=sum(phase=="TEST"); report.nHabituation=sum(phase=="HABITUATION");
report.nFLS=sum(typeShort=="FLS"); report.nFLD=sum(typeShort=="FLD");
report.nRLS=sum(typeShort=="RLS"); report.nRLD=sum(typeShort=="RLD");
report.nStimPaired=sum(stimFlag); report.nStimGS=sum(isMatchedGS); report.nStimGD=sum(isStimGD);
report.nStimHab=sum(stimFlag & phase=="HABITUATION"); report.nPostStimBuffers=sum(postStimBuffer);
report.bufferRuleDetected=bufferRuleDetected; report.nBufferRuleViolations=sum(bufferRuleViolation(stimIdx));
report.nExplicitBufferMessages=sum(explicitBufferMessage);

out=struct();
out.tone=struct('t_msg_ms',tmsg,'t_onset_ms',tmsg,'code',code,'physical',phys,'nTrials',nTrials);
out.trialInfo=trialInfo; out.report=report;

if verbose
    fprintf('\n[EyeLink ER-taVNS LG]\n');
    fprintf('  Complete 5-tone trials: %d / %d\n',report.completeToneTrials,expectedTrials);
    fprintf('  T5=%d | TEST=%d | HAB=%d\n',report.nT5,report.nTest,report.nHabituation);
    fprintf('  FLS=%d FLD=%d RLS=%d RLD=%d\n',report.nFLS,report.nFLD,report.nRLS,report.nRLD);
    fprintf('  paired stim=%d | matched GS=%d | stimulated GD=%d | HAB stim=%d | derived buffers=%d\n', ...
        report.nStimPaired,report.nStimGS,report.nStimGD,report.nStimHab,report.nPostStimBuffers);
    fprintf('  N+1 buffer rule detected: %d | violations=%d | historical explicit buffer messages=%d\n', ...
        report.bufferRuleDetected,report.nBufferRuleViolations,report.nExplicitBufferMessages);
    fprintf('  Duplicate tone slots: %d\n',report.duplicateToneSlots);
    if report.completeToneTrials~=expectedTrials, warning('[EyeLink LG] Not all expected trials have five tone messages.'); end
    if report.duplicateToneSlots>0, warning('[EyeLink LG] Duplicate trial/tone message slots detected.'); end
    if any([report.nFLS report.nFLD report.nRLS report.nRLD]~=[169 169 36 36])
        warning('[EyeLink LG] Condition counts differ from 169/169/36/36.');
    end
    if report.nStimPaired~=72 || report.nStimGS~=36 || report.nStimGD~=36 || report.nStimHab~=0
        warning('[EyeLink LG] Stimulation counts differ from expected 72 stim = 36 GS + 36 GD, 0 HAB stim.');
    end
    if report.bufferRuleDetected && report.nPostStimBuffers~=72
        warning('[EyeLink LG] N+1 buffer rule was detected but the derived buffer count is not 72.');
    elseif ~report.bufferRuleDetected && report.nStimPaired>0
        warning('[EyeLink LG] Mandatory N+1 buffer rule not detected; automatic buffer exclusion is disabled for this recording.');
    end
end
end
