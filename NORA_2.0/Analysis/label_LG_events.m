function [EEG, report] = label_LG_events(EEG, varargin)
% LABEL_LG_EVENTS
% Decode the CURRENT ER-taVNS 4-bit BrainVision marker stream.
%
% Current acquisition markers:
%   S1 = FLS condition label
%   S2 = FLD condition label
%   S3 = RLS condition label
%   S4 = RLD condition label
%   S5 = T1
%   S6 = T2
%   S7 = T3
%   S8 = T4
%   S9 = T5 (PRIMARY EEG epoch / stimulation anchor)
%  S10 = stimulation-paired trial label (GD OR matched GS)
%  S11..S14 = E2 visual markers; S15 = E2 response
%
% IMPORTANT:
%   S1..S4 and S10 are emitted a few milliseconds AFTER S9.  Therefore the
%   T5 marker S9 must remain the epoch anchor.  This function reads the
%   following S1..S4 marker, relabels the ORIGINAL S9 event as FLS/FLD/RLS/RLD,
%   and annotates whether an S10 followed it.  No event latencies are shifted.
%
% The first 25 auditory patterns in each 135-pattern block are habituation.
% Habituation uses the same FLS/FLD condition labels but is never stimulated.
%
% Added fields on the auditory tone events include:
%   lg_code, lg_tone_position, lg_trial, lg_block, lg_trial_in_block
%   lg_phase, lg_type_short, lg_type, lg_local_short/lg_local
%   lg_global_short/lg_global, lg_is_anchor
%   lg_stim_paired, lg_stim_kind, lg_post_stim_buffer, lg_exclude_primary
%   lg_condition_marker_latency, lg_stim_marker_latency
%
% The condition-label event itself is renamed COND_FLS/COND_FLD/COND_RLS/
% COND_RLD.  S10 is renamed STIM_PAIRED.  E2 S11..S15 are untouched.

p = inputParser;
p.addParameter('verbose', true, @(x)islogical(x) && isscalar(x));
p.addParameter('pairWindow_ms', 50, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
verbose = p.Results.verbose;
pairWindow_ms = p.Results.pairWindow_ms;

if ~isfield(EEG,'event') || isempty(EEG.event)
    error('EEG.event is empty.');
end

FS = double(EEG.srate);
pairWindowSamples = pairWindow_ms * FS / 1000;

report = struct();
report.nT1=0; report.nT2=0; report.nT3=0; report.nT4=0; report.nT5=0;
report.nFLS=0; report.nFLD=0; report.nRLS=0; report.nRLD=0;
report.nConditionMarkers=0; report.nStimPaired=0;
report.nStimGS=0; report.nStimGD=0; report.nStimHabituation=0;
report.nMissingCondition=0; report.nCompleteFiveToneTrials=0;
report.nHabAnchors=0; report.nTestAnchors=0;
report.nConsecutiveStimPairs=0; report.nPostStimBuffers=0;
report.conditionOffset_ms=[]; report.stimOffset_ms=[];

% Preserve original marker types and decode numeric codes once.
nEv = numel(EEG.event);
rawCode = nan(1,nEv);
for i=1:nEv
    if ~isfield(EEG.event(i),'type'), continue; end
    if ~isfield(EEG.event(i),'original_type') || isempty(EEG.event(i).original_type)
        EEG.event(i).original_type = EEG.event(i).type;
    end
    rawCode(i)=parseStimCode(EEG.event(i).type);
end

latAll = double([EEG.event.latency]);
[~,orderAll] = sort(latAll);

% Rename the simple raw markers first.  S9 is temporary T5 until its
% following condition label has been decoded.
for i=1:nEv
    c=rawCode(i);
    switch c
        case 5, EEG.event(i).type='T1'; report.nT1=report.nT1+1;
        case 6, EEG.event(i).type='T2'; report.nT2=report.nT2+1;
        case 7, EEG.event(i).type='T3'; report.nT3=report.nT3+1;
        case 8, EEG.event(i).type='T4'; report.nT4=report.nT4+1;
        case 9, EEG.event(i).type='T5'; report.nT5=report.nT5+1;
        case 1, EEG.event(i).type='COND_FLS'; report.nConditionMarkers=report.nConditionMarkers+1;
        case 2, EEG.event(i).type='COND_FLD'; report.nConditionMarkers=report.nConditionMarkers+1;
        case 3, EEG.event(i).type='COND_RLS'; report.nConditionMarkers=report.nConditionMarkers+1;
        case 4, EEG.event(i).type='COND_RLD'; report.nConditionMarkers=report.nConditionMarkers+1;
        case 10, EEG.event(i).type='STIM_PAIRED'; report.nStimPaired=report.nStimPaired+1;
    end
end

trialToneIdx = nan(1,4);
trialCounter = 0;
stimByTrial = false(1,0); bufferByTrial = false(1,0);
globalByTrial = strings(1,0);

for oo=1:numel(orderAll)
    idx=orderAll(oo);
    c=rawCode(idx);

    if c==5
        trialToneIdx(:)=NaN; trialToneIdx(1)=idx;
    elseif c==6
        if isfinite(trialToneIdx(1)), trialToneIdx(2)=idx; end
    elseif c==7
        if all(isfinite(trialToneIdx(1:2))), trialToneIdx(3)=idx; end
    elseif c==8
        if all(isfinite(trialToneIdx(1:3))), trialToneIdx(4)=idx; end
    elseif c==9
        trialCounter=trialCounter+1;
        t5Lat=double(EEG.event(idx).latency);

        % Search only forward from S9, within 50 ms, and stop before the
        % next S5/T1 if one appears unusually early.
        condIdx=[]; stimIdx=[];
        for pp=(oo+1):numel(orderAll)
            j=orderAll(pp);
            dt=double(EEG.event(j).latency)-t5Lat;
            if dt<0, continue; end
            if dt>pairWindowSamples, break; end
            cj=rawCode(j);
            if cj==5, break; end
            if ismember(cj,1:4)
                if isempty(condIdx), condIdx=j; end
            end
            if isempty(stimIdx) && cj==10, stimIdx=j; end
        end

        if isempty(condIdx)
            report.nMissingCondition=report.nMissingCondition+1;
            trialToneIdx(:)=NaN;
            continue;
        end

        condCode=rawCode(condIdx);
        % Buffer status is NOT encoded by an extra TTL. It is assigned after
        % all trials are decoded from the deterministic rule: N+1 after S10.
        isPostStimBuffer=false;
        [short,longType,localShort,localLong,globalShort,globalLong]=conditionInfo(condCode);
        EEG.event(idx).type=short; % S9 remains the latency anchor, now condition-labeled.

        inBlock=mod(trialCounter-1,205)+1;
        blockNo=ceil(trialCounter/205);
        if inBlock<=25
            phase='HABITUATION'; report.nHabAnchors=report.nHabAnchors+1;
        else
            phase='TEST'; report.nTestAnchors=report.nTestAnchors+1;
        end
        stimPaired=~isempty(stimIdx);
        if stimPaired && strcmp(globalShort,'GS')
            stimKind='MATCHED_GS'; report.nStimGS=report.nStimGS+1;
        elseif stimPaired && strcmp(globalShort,'GD')
            stimKind='STIMULATED_GD'; report.nStimGD=report.nStimGD+1;
        else
            stimKind='NONE';
        end
        if stimPaired && strcmp(phase,'HABITUATION')
            report.nStimHabituation=report.nStimHabituation+1;
        end
        stimByTrial(trialCounter)=stimPaired; %#ok<AGROW>
        bufferByTrial(trialCounter)=false; %#ok<AGROW>
        globalByTrial(trialCounter)=string(globalShort); %#ok<AGROW>

        condOffset_ms=1000*(double(EEG.event(condIdx).latency)-t5Lat)/FS;
        report.conditionOffset_ms(end+1)=condOffset_ms; %#ok<AGROW>
        if stimPaired
            stimOffset_ms=1000*(double(EEG.event(stimIdx).latency)-t5Lat)/FS;
            report.stimOffset_ms(end+1)=stimOffset_ms; %#ok<AGROW>
        else
            stimOffset_ms=NaN;
        end

        % Metadata applied to T5 and the preceding T1-T4 events.
        toneIdx=[trialToneIdx idx];
        if all(isfinite(trialToneIdx)), report.nCompleteFiveToneTrials=report.nCompleteFiveToneTrials+1; end
        for pos=1:numel(toneIdx)
            j=toneIdx(pos);
            if ~isfinite(j), continue; end
            EEG.event(j).lg_code=rawCode(j);
            EEG.event(j).lg_tone_position=pos;
            EEG.event(j).lg_trial=trialCounter;
            EEG.event(j).lg_block=blockNo;
            EEG.event(j).lg_trial_in_block=inBlock;
            EEG.event(j).lg_phase=phase;
            EEG.event(j).lg_type_short=short;
            EEG.event(j).lg_type=longType;
            EEG.event(j).lg_local_short=localShort;
            EEG.event(j).lg_local=localLong;
            EEG.event(j).lg_global_short=globalShort;
            EEG.event(j).lg_global=globalLong;
            EEG.event(j).lg_is_anchor=(pos==5);
            EEG.event(j).lg_stim_paired=stimPaired;
            EEG.event(j).lg_stim_kind=stimKind;
            EEG.event(j).lg_post_stim_buffer=isPostStimBuffer;
            EEG.event(j).lg_exclude_primary=isPostStimBuffer;
            EEG.event(j).lg_condition_marker_latency=double(EEG.event(condIdx).latency);
            EEG.event(j).lg_stim_marker_latency=NaN;
            if stimPaired, EEG.event(j).lg_stim_marker_latency=double(EEG.event(stimIdx).latency); end
            EEG.event(j).lg_condition_offset_ms=condOffset_ms;
            EEG.event(j).lg_stim_offset_ms=stimOffset_ms;
        end

        % Also annotate the descriptive condition/S10 events themselves.
        for j=[condIdx stimIdx]
            if isempty(j), continue; end
            EEG.event(j).lg_trial=trialCounter;
            EEG.event(j).lg_block=blockNo;
            EEG.event(j).lg_trial_in_block=inBlock;
            EEG.event(j).lg_phase=phase;
            EEG.event(j).lg_type_short=short;
            EEG.event(j).lg_type=longType;
            EEG.event(j).lg_local_short=localShort;
            EEG.event(j).lg_local=localLong;
            EEG.event(j).lg_global_short=globalShort;
            EEG.event(j).lg_global=globalLong;
            EEG.event(j).lg_stim_paired=stimPaired;
            EEG.event(j).lg_stim_kind=stimKind;
            EEG.event(j).lg_post_stim_buffer=isPostStimBuffer;
            EEG.event(j).lg_exclude_primary=isPostStimBuffer;
            EEG.event(j).lg_is_anchor=false;
        end


        switch short
            case 'FLS', report.nFLS=report.nFLS+1;
            case 'FLD', report.nFLD=report.nFLD+1;
            case 'RLS', report.nRLS=report.nRLS+1;
            case 'RLD', report.nRLD=report.nRLD+1;
        end

        trialToneIdx(:)=NaN;
    end
end

if numel(stimByTrial)>1
    report.nConsecutiveStimPairs=sum(stimByTrial(1:end-1)&stimByTrial(2:end));
end
report.nTrials=trialCounter;

% ---------------------------------------------------------------------
% Robust buffer derivation: the acquisition rule is that every S10 trial N
% is followed immediately by an UNSTIMULATED GS trial N+1. No dedicated
% buffer TTL is required. First verify the complete session obeys that rule.
% This also keeps legacy pre-buffer pilots safe: if the rule is not obeyed
% globally, no N+1 exclusions are inferred from S10 alone.
report.bufferRuleDetected=false;
report.nStimRuleViolations=0;
if ~isempty(stimByTrial) && any(stimByTrial)
    stimTrials=find(stimByTrial);
    bad=false(size(stimTrials));
    for q=1:numel(stimTrials)
        n=stimTrials(q);
        bad(q) = n>=numel(stimByTrial) || stimByTrial(n+1) || globalByTrial(n+1)~="GS";
    end
    report.nStimRuleViolations=sum(bad);
    report.bufferRuleDetected=~any(bad);
    if report.bufferRuleDetected
        bufferByTrial(2:end)=stimByTrial(1:end-1);
    end
end
report.nPostStimBuffers=sum(bufferByTrial);
report.nStimWithoutFollowingBuffer=report.nStimRuleViolations;

% Apply derived buffer metadata to every event belonging to each N+1 trial.
% The T5 anchor then carries this flag into EEG.epoch via set_LG_epoch_labels.
for tr=find(bufferByTrial)
    ev=find(arrayfun(@(x)isfield(x,'lg_trial') && isequal(x.lg_trial,tr),EEG.event));
    for q=ev
        EEG.event(q).lg_post_stim_buffer=true;
        EEG.event(q).lg_exclude_primary=true;
    end
end

EEG=eeg_checkset(EEG,'eventconsistency');
EEG.etc.LG_trigger_report=report;

if verbose
    fprintf('\n[ER-taVNS LG trigger/label check]\n');
    fprintf('  T1=%d T2=%d T3=%d T4=%d T5=%d\n',report.nT1,report.nT2,report.nT3,report.nT4,report.nT5);
    fprintf('  FLS=%d FLD=%d RLS=%d RLD=%d\n',report.nFLS,report.nFLD,report.nRLS,report.nRLD);
    fprintf('  S10 paired=%d | matched GS=%d | stimulated GD=%d | HAB stim=%d\n', ...
        report.nStimPaired,report.nStimGS,report.nStimGD,report.nStimHabituation);
    fprintf('  N+1 buffers derived=%d | rule violations=%d | rule detected=%d\n', ...
        report.nPostStimBuffers,report.nStimRuleViolations,report.bufferRuleDetected);
    fprintf('  Complete 5-tone chains=%d | missing condition after S9=%d\n', ...
        report.nCompleteFiveToneTrials,report.nMissingCondition);
    if ~isempty(report.conditionOffset_ms)
        fprintf('  Median S9->condition offset = %.2f ms\n',median(report.conditionOffset_ms,'omitnan'));
    end
    if ~isempty(report.stimOffset_ms)
        fprintf('  Median S9->S10 offset      = %.2f ms\n',median(report.stimOffset_ms,'omitnan'));
    end

    if any([report.nT1 report.nT2 report.nT3 report.nT4 report.nT5]~=410)
        warning('[LG] Expected 410 occurrences of each T1..T5 marker.');
    end
    if any([report.nFLS report.nFLD report.nRLS report.nRLD]~=[169 169 36 36])
        warning('[LG] Expected FLS/FLD/RLS/RLD = 169/169/36/36.');
    end
    if report.nStimPaired~=72 || report.nStimGS~=36 || report.nStimGD~=36
        warning('[LG] Expected S10 = 72 total (36 GS + 36 GD).');
    end
    if report.nStimHabituation~=0
        warning('[LG] Stimulation marker found during habituation. Expected 0.');
    end
    if report.nConsecutiveStimPairs~=0
        warning('[LG] Consecutive stimulation-paired trials detected.');
    end
    if report.nStimPaired>0 && (~report.bufferRuleDetected || report.nPostStimBuffers~=72)
        warning('[LG] N+1 buffer rule was not detected cleanly. Legacy pre-buffer data are not auto-excluded; current-design data should have 72 derived buffers.');
    end
end
end

% =====================================================================
function code=parseStimCode(x)
code=NaN;
if isnumeric(x)&&isscalar(x)&&isfinite(x), code=double(x); return; end
if iscell(x), if isempty(x),return;end, x=x{1}; end
s=strtrim(char(string(x)));
tok=regexp(s,'(?i)(?:^|/|\s)S\s*([0-9]+)\s*$','tokens','once');
if ~isempty(tok), code=str2double(tok{1}); return; end
v=str2double(s); if isfinite(v), code=v; end
end

% =====================================================================
function [short,longType,localShort,localLong,globalShort,globalLong]=conditionInfo(code)
switch code
    case 1
        short='FLS'; longType='Frequent Local Standard';
        localShort='LS'; localLong='Local Standard'; globalShort='GS'; globalLong='Global Standard';
    case 2
        short='FLD'; longType='Frequent Local Deviant';
        localShort='LD'; localLong='Local Deviant'; globalShort='GS'; globalLong='Global Standard';
    case 3
        short='RLS'; longType='Rare Local Standard';
        localShort='LS'; localLong='Local Standard'; globalShort='GD'; globalLong='Global Deviant';
    case 4
        short='RLD'; longType='Rare Local Deviant';
        localShort='LD'; localLong='Local Deviant'; globalShort='GD'; globalLong='Global Deviant';
    otherwise
        error('Unknown LG condition code %g.',code);
end
end
