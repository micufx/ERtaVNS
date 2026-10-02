function results = preprocess_pupil_LG(lines, lg, fs, dt, pad_ms, relTimes_out, b_lp, a_lp, blinkLimit, varargin)
% PREPROCESS_PUPIL_LG
% Corrected pupil preprocessing for the new Exp_01 / Exp_02 marker scheme.
% Keeps your original pipeline:
% - blink mask
% - fillmissing linear across blinks
% - low-pass filter
% - epoch around TONE5 (with interp1 to uniform time axis)
% - baseline correction (relative to TONE5 or TONE1)
% - blink/NaN exclusion
% - z-transform
%
% OUTPUT epoch is ALWAYS relTimes_out (e.g., -500..2000)

p = inputParser;
p.addParameter('baselineAnchor', 'TONE5', @(s)ischar(s)||isstring(s));
p.addParameter('baselineWin_ms', [-500 0], @(x)isnumeric(x)&&numel(x)==2);
p.parse(varargin{:});
baselineAnchor = upper(string(p.Results.baselineAnchor));
baselineWin_ms = p.Results.baselineWin_ms(:)';

if ~isstruct(lg) || ~isfield(lg,'tone') || ~isfield(lg,'trialInfo')
    error('lg must be the output of parse_LG_eyelink.');
end
tone = lg.tone;
trialInfo = lg.trialInfo;
t = tone.t_onset_ms;

% --- Decide internal big window (only to ensure baseline exists) ---
epochWin_out_ms = [relTimes_out(1) relTimes_out(end)];

needExtraPre = 0;
if baselineAnchor=="TONE1"
    have15 = isfinite(t(:,1)) & isfinite(t(:,5));
    shift15 = t(have15,1) - t(have15,5);
    if ~isempty(shift15)
        baseStartMin = min(shift15 + baselineWin_ms(1));
        if baseStartMin < epochWin_out_ms(1)
            needExtraPre = (epochWin_out_ms(1) - baseStartMin) + 50; % safety 50 ms
        end
    end
end

% Make big window aligned to dt grid (IMPORTANT!)
pre_ms_big  = epochWin_out_ms(1) - needExtraPre - 100;
post_ms_big = epochWin_out_ms(2) + 200;

pre_ms_big  = dt * floor(pre_ms_big/dt);
post_ms_big = dt * ceil(post_ms_big/dt);

relTimes_big = pre_ms_big : dt : post_ms_big;

% Compute index mapping from output times to big times (NO ismember)
idx_out = round((relTimes_out - pre_ms_big)/dt) + 1;
if any(idx_out < 1) || any(idx_out > numel(relTimes_big))
    error('Output epoch indices fall outside internal big epoch. Check dt/grid.');
end

%% --- Parse sample lines ---
sampleTimes=[]; sampleX=[]; sampleY=[]; samplePupil=[]; sampleFlag=[];
for k=1:numel(lines)
    L=strtrim(lines{k});
    if isempty(L), continue; end
    if isstrprop(L(1),'digit')
        vals = sscanf(L,'%f',5);
        if numel(vals)>=5
            sampleTimes(end+1,1)=vals(1);
            sampleX(end+1,1)=vals(2);
            sampleY(end+1,1)=vals(3);
            samplePupil(end+1,1)=vals(4);
            sampleFlag(end+1,1)=vals(5);
        end
    end
end
if isempty(sampleTimes), error('No sample lines parsed.'); end

%% --- Blink intervals ---
blinkStart=[]; blinkEnd=[];
for k=1:numel(lines)
    L=strtrim(lines{k});
    if isempty(L), continue; end
    tok = regexp(L,'^EBLINK\s+\w+\s+(\d+)\s+(\d+)\s+(\d+)','tokens','once');
    if ~isempty(tok)
        blinkStart(end+1,1)=str2double(tok{1});
        blinkEnd(end+1,1)=str2double(tok{2});
    end
end

%% --- Continuous preprocessing (same as your original) ---
blinkMask_cont = false(size(sampleTimes));
for b=1:numel(blinkStart)
    bStartPad = blinkStart(b) - pad_ms;
    bEndPad   = blinkEnd(b) + pad_ms;
    blinkMask_cont(sampleTimes>=bStartPad & sampleTimes<=bEndPad)=true;
end

pupil_cont = samplePupil;
pupil_cont(blinkMask_cont)=NaN;
pupil_cont = fillmissing(pupil_cont,'linear');
pupil_filt = filtfilt(b_lp,a_lp,pupil_cont);

%% --- Epoch around TONE5 using parsed LG trial metadata ---
% The corrected acquisition writes all labels at T5 onset. We therefore
% do not search for a later AUDIO_TRIGGER and do not apply a 50-ms shift.

type=strings(0,1); effect=strings(0,1); sound=strings(0,1);
patternStr=strings(0,1); modality=strings(0,1);
timeVal=zeros(0,1); block=zeros(0,1); trialInFile=zeros(0,1);

X_GAZE_epoch=cell(0,1); Y_GAZE_epoch=cell(0,1);
PUPIL_epoch=cell(0,1); EVENTFLAG_epoch=cell(0,1); BLINK_epoch=cell(0,1);
stimPaired=false(0,1); postStimBuffer=false(0,1); stimKind=strings(0,1); sessionCondition=strings(0,1); visitNum=zeros(0,1);

skippedWindow=0;
skippedNoLabel=0;
skippedNoSamples=0;

skippedNoSamples_sound=strings(0,1);
skippedNoSamples_block=zeros(0,1);
skippedNoSamples_time=zeros(0,1);

skippedWindow_sound=strings(0,1);
skippedWindow_block=zeros(0,1);
skippedWindow_time=zeros(0,1);

nTrialsAvailable=min(size(t,1),height(trialInfo));

for trialIdx=1:nTrialsAvailable

    if ~isfinite(t(trialIdx,5))
        skippedNoLabel=skippedNoLabel+1;
        continue;
    end

    t5=t(trialIdx,5);

    rawType=upper(string(trialInfo.raw_type_short(trialIdx)));
    normType=upper(string(trialInfo.type_short(trialIdx)));
    localShort=upper(string(trialInfo.local_short(trialIdx)));
    globalShort=upper(string(trialInfo.global_short(trialIdx)));
    phaseStr=upper(string(trialInfo.phase(trialIdx)));

    thisStim=false; thisStimKind="NONE"; thisSession=""; thisVisit=NaN;
    if ismember('stim_paired',trialInfo.Properties.VariableNames), thisStim=logical(trialInfo.stim_paired(trialIdx)); end
    thisPostStimBuffer=false;
    if ismember('is_post_stim_buffer',trialInfo.Properties.VariableNames), thisPostStimBuffer=logical(trialInfo.is_post_stim_buffer(trialIdx)); end
    if ismember('is_matched_gs_control',trialInfo.Properties.VariableNames) && trialInfo.is_matched_gs_control(trialIdx)
        thisStimKind="MATCHED_GS";
    elseif ismember('is_stimulated_gd',trialInfo.Properties.VariableNames) && trialInfo.is_stimulated_gd(trialIdx)
        thisStimKind="STIMULATED_GD";
    end
    if ismember('session_condition',trialInfo.Properties.VariableNames), thisSession=string(trialInfo.session_condition(trialIdx)); end
    if ismember('visit',trialInfo.Properties.VariableNames), thisVisit=trialInfo.visit(trialIdx); end

    if normType=="FLS"
        typeStr="Frequent Local Standard";
    elseif normType=="FLD"
        typeStr="Frequent Local Deviant";
    elseif normType=="RLS"
        typeStr="Rare Local Standard";
    elseif normType=="RLD"
        typeStr="Rare Local Deviant";
    else
        skippedNoLabel=skippedNoLabel+1;
        continue;
    end

    if localShort=="LS"
        soundStr="Local Standard";
        patt="[1,1,1,1,1]";
    elseif localShort=="LD"
        soundStr="Local Deviant";
        patt="[1,1,1,1,2]";
    else
        skippedNoLabel=skippedNoLabel+1;
        continue;
    end

    if globalShort=="GS"
        effectStr="Global Standard";
    elseif globalShort=="GD"
        effectStr="Global Deviant";
    else
        skippedNoLabel=skippedNoLabel+1;
        continue;
    end

    blockId=trialInfo.block(trialIdx);
    if phaseStr=="HABITUATION"
        modalityStr="habituation";
    else
        modalityStr="test";
    end

    % Window in absolute EyeLink time.
    wStart=t5+relTimes_big(1);
    wEnd=t5+relTimes_big(end);

    if wStart<sampleTimes(1) || wEnd>sampleTimes(end)
        skippedWindow=skippedWindow+1;
        skippedWindow_sound(end+1,1)=soundStr;
        skippedWindow_block(end+1,1)=blockId;
        skippedWindow_time(end+1,1)=t5;
        continue;
    end

    inWin=(sampleTimes>=wStart)&(sampleTimes<=wEnd);

    if nnz(inWin)<2
        skippedNoSamples=skippedNoSamples+1;
        skippedNoSamples_sound(end+1,1)=soundStr;
        skippedNoSamples_block(end+1,1)=blockId;
        skippedNoSamples_time(end+1,1)=t5;
        continue;
    end

    tsLocal=sampleTimes(inWin)-t5;
    xLocal=sampleX(inWin);
    yLocal=sampleY(inWin);
    pupLocal=pupil_filt(inWin);
    flgLocal=sampleFlag(inWin);
    blnLocal=blinkMask_cont(inWin);

    goodT=isfinite(tsLocal);
    tsLocal=tsLocal(goodT);
    xLocal=xLocal(goodT);
    yLocal=yLocal(goodT);
    pupLocal=pupLocal(goodT);
    flgLocal=flgLocal(goodT);
    blnLocal=blnLocal(goodT);

    if numel(tsLocal)<2
        skippedNoSamples=skippedNoSamples+1;
        skippedNoSamples_sound(end+1,1)=soundStr;
        skippedNoSamples_block(end+1,1)=blockId;
        skippedNoSamples_time(end+1,1)=t5;
        continue;
    end

    [tsS,sIdx]=sort(tsLocal);
    xLocal=xLocal(sIdx);
    yLocal=yLocal(sIdx);
    pupLocal=pupLocal(sIdx);
    flgLocal=flgLocal(sIdx);
    blnLocal=blnLocal(sIdx);

    [tsLocal,uIdx]=unique(tsS,'stable');
    xLocal=xLocal(uIdx);
    yLocal=yLocal(uIdx);
    pupLocal=pupLocal(uIdx);
    flgLocal=flgLocal(uIdx);
    blnLocal=blnLocal(uIdx);

    if numel(tsLocal)<2
        skippedNoSamples=skippedNoSamples+1;
        skippedNoSamples_sound(end+1,1)=soundStr;
        skippedNoSamples_block(end+1,1)=blockId;
        skippedNoSamples_time(end+1,1)=t5;
        continue;
    end

    xBig=interp1(tsLocal,xLocal,relTimes_big,'linear',NaN);
    yBig=interp1(tsLocal,yLocal,relTimes_big,'linear',NaN);
    pupBig=interp1(tsLocal,pupLocal,relTimes_big,'linear',NaN);
    flgBig=interp1(tsLocal,flgLocal,relTimes_big,'nearest',NaN);
    blnBig=interp1(tsLocal,double(blnLocal),relTimes_big,'nearest',0)>0;

    % Baseline correction on the BIG axis: unchanged from the old pipeline.
    if baselineAnchor=="TONE5"
        wRel=baselineWin_ms;
    else
        if ~isfinite(t(trialIdx,1)) || ~isfinite(t(trialIdx,5))
            wRel=[NaN NaN];
        else
            wRel=(t(trialIdx,1)-t(trialIdx,5))+baselineWin_ms;
        end
    end

    if all(isfinite(wRel))
        idxB=(relTimes_big>=wRel(1))&(relTimes_big<=wRel(2));
        if nnz(idxB)>=5
            b=mean(pupBig(idxB),'omitnan');
            pupBig=pupBig-b;
        end
    end

    xEpoch=xBig(idx_out);
    yEpoch=yBig(idx_out);
    pupEpoch=pupBig(idx_out);
    flgEpoch=flgBig(idx_out);
    blnEpoch=blnBig(idx_out);

    type(end+1,1)=typeStr;
    effect(end+1,1)=effectStr;
    sound(end+1,1)=soundStr;
    patternStr(end+1,1)=patt;
    modality(end+1,1)=modalityStr;
    timeVal(end+1,1)=t5;
    block(end+1,1)=blockId;
    trialInFile(end+1,1)=trialIdx;
    stimPaired(end+1,1)=thisStim;
    postStimBuffer(end+1,1)=thisPostStimBuffer;
    stimKind(end+1,1)=thisStimKind;
    sessionCondition(end+1,1)=thisSession;
    visitNum(end+1,1)=thisVisit;

    X_GAZE_epoch{end+1,1}=xEpoch;
    Y_GAZE_epoch{end+1,1}=yEpoch;
    PUPIL_epoch{end+1,1}=pupEpoch;
    EVENTFLAG_epoch{end+1,1}=flgEpoch;
    BLINK_epoch{end+1,1}=blnEpoch;
end

if skippedNoLabel>0
    warning('[PUPIL] %d trials lacked complete LG T5 labels and were skipped.',skippedNoLabel);
end

dataTable = table(type,effect,sound,patternStr,modality,timeVal,block,trialInFile, ...
    stimPaired,postStimBuffer,stimKind,sessionCondition,visitNum, ...
    X_GAZE_epoch,Y_GAZE_epoch,PUPIL_epoch,EVENTFLAG_epoch,BLINK_epoch, ...
    'VariableNames', {'type','effect','sound','pattern','modality','time','block','trialInFile', ...
                      'stim_paired','is_post_stim_buffer','stim_kind','session_condition','visit', ...
                      'X_GAZE','Y_GAZE','PUPIL','EVENT_FLAG','BLINK'});

% ---- Exclusions (same logic you already had) ----
if height(dataTable)==0
    results.dataTableClean = dataTable;
    results.removedTrialsSummaryAll = table();
    results.relTimes_out = relTimes_out;
    return;
end

nSamples = numel(dataTable.PUPIL{1});
keepMask = true(height(dataTable),1);
removedPostStimBuffer = find(dataTable.is_post_stim_buffer);
keepMask(removedPostStimBuffer)=false;
blinkOK = true(height(dataTable),1);
for tr=1:height(dataTable)
    blinkOK(tr) = (sum(dataTable.BLINK{tr})/nSamples) <= blinkLimit;
end
removedBlink = find(~blinkOK & ~dataTable.is_post_stim_buffer);
keepMask = keepMask & blinkOK;

hasNaN = cellfun(@(x) any(isnan(x)), dataTable.PUPIL);
removedNaN = find(hasNaN & keepMask);

finalKeepMask = keepMask;
finalKeepMask(removedNaN)=false;

RemovalReason = strings(height(dataTable),1);
RemovalReason(removedPostStimBuffer)="Post-stimulation GS buffer (predefined exclusion)";
RemovalReason(removedBlink)="Blink >20%";
RemovalReason(removedNaN)="Trial incomplete (missing data)";

removedTrialsSummary = table( ...
    find(~finalKeepMask), ...
    RemovalReason(~finalKeepMask), ...
    dataTable.sound(~finalKeepMask), ...
    dataTable.block(~finalKeepMask), ...
    dataTable.time(~finalKeepMask), ...
    'VariableNames', {'TrialIndex','Reason','SoundType','Block','ToneTime'} );

if skippedWindow>0
    removedWindowSummary = table( ...
        nan(skippedWindow,1), ...
        repmat("Not fully inside recording window",skippedWindow,1), ...
        skippedWindow_sound, skippedWindow_block, skippedWindow_time, ...
        'VariableNames', {'TrialIndex','Reason','SoundType','Block','ToneTime'} );
    removedTrialsSummaryAll = [removedTrialsSummary; removedWindowSummary];
else
    removedTrialsSummaryAll = removedTrialsSummary;
end

% NEW (minimal): add "No samples in epoch window" skips into the SAME summary table
if skippedNoSamples > 0
    removedNoSamplesSummary = table( ...
        nan(skippedNoSamples,1), ...
        repmat("No samples in epoch window",skippedNoSamples,1), ...
        skippedNoSamples_sound, skippedNoSamples_block, skippedNoSamples_time, ...
        'VariableNames', {'TrialIndex','Reason','SoundType','Block','ToneTime'} );
    removedTrialsSummaryAll = [removedTrialsSummaryAll; removedNoSamplesSummary];
end

dataTableClean = dataTable(finalKeepMask,:);

% ---- Z-transform ----
nTrials = height(dataTableClean);
if nTrials > 0
    pupilMat = nan(nTrials, numel(relTimes_out));
    for tr=1:nTrials
        pupilMat(tr,:) = dataTableClean.PUPIL{tr}(:).';
    end
    globalMean = mean(pupilMat(:),'omitnan');
    globalSD   = std(pupilMat(:),'omitnan');
    pupilZmat  = (pupilMat-globalMean)/globalSD;

    PUPIL_Z = cell(nTrials,1);
    for tr=1:nTrials
        PUPIL_Z{tr} = pupilZmat(tr,:);
    end
    dataTableClean.PUPIL_Z = PUPIL_Z;
else
    dataTableClean.PUPIL_Z = cell(0,1);
end

results.dataTableClean          = dataTableClean;
results.removedTrialsSummaryAll = removedTrialsSummaryAll;
results.relTimes_out            = relTimes_out;

end

function out = iff(cond,a,b)
if cond, out=a; else, out=b; end
end
