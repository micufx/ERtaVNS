clc; clearvars; close all;

%% CONTROL PHASE 
% Keep this experiment folder first on the MATLAB path (io64 files live here).
THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(THIS_SCRIPT_DIR)
    addpath(THIS_SCRIPT_DIR, '-begin');
end
try, Eyelink('Shutdown'); catch, end
clear mex
P = ERtaVNS_parameters(); % central editable stimulation/timing configuration
% Complete NI packet is prequeued before each burst; no DataRequired refill callback.

%% ================================================================
%  VNS CONTROL PHASE | 2-SESSION STIMULATION-ONLY CONTROL
%
%  Per visit:
%    - 30 bursts of the visit condition ONLY
%      ACTIVE visit: left cymba conchae
%      SHAM visit:   left earlobe
%    - burst = 1500 ms, 25 Hz, 250 us
%    - OFF interval = 5.0-5.3 s, continuous random jitter
%    - fixation only; no auditory tones
%    - one DS5 / NI Dev2 AO0 pathway
%
%  UNIFIED FILE-NAMING CONVENTION:
%      Visit 1 control: S01V1C.edf / S01V1C.mat
%      Visit 2 control: S01V2C.edf / S01V2C.mat
%  The SAME base stem is used for EyeLink, local MATLAB data and the
%  manually named BrainVision recording. Visit folders are retained only
%  for organization; the filename itself already encodes the visit.
%% ================================================================

%% ---------------- PARTICIPANT + VISIT INFORMATION ----------------
prompt={'Participant number / ID (e.g., S01):','Visit number (1 or 2):'};
dlg_title='ER-taVNS | stimulation-only control';
def={'S01','1'};
answer=inputdlg(prompt,dlg_title,1,def);
if isempty(answer), error('Control experiment cancelled before start.'); end
participantInput=upper(strtrim(answer{1}));
tok=regexp(participantInput,'^S?(\d+)$','tokens','once');
if isempty(tok), error('Participant ID must have the form S01 or a positive integer.'); end
participantNumber=str2double(tok{1});
visitNumber=str2double(strtrim(answer{2}));
if ~isfinite(participantNumber) || participantNumber<1 || participantNumber~=round(participantNumber)
    error('Participant number must be a positive integer.');
end
if ~ismember(visitNumber,[1 2]), error('Visit number must be 1 or 2.'); end
[participantNumber,participantID,visitNumber]=ERtaVNS_resolve_participant(participantNumber,visitNumber,'C',THIS_SCRIPT_DIR);
assignment=ERtaVNS_get_assignment(participantNumber,visitNumber);
sessionCondition=assignment.sessionCondition;
sessionSite=ERtaVNS_condition_site(sessionCondition);
% Unified EyeLink/MATLAB/BrainVision base name (<= 8 EyeLink characters).
baseFileStem=sprintf('%sV%dC',participantID,visitNumber);
if length(baseFileStem)>8, error('Generated base file name exceeds the EyeLink 8-character limit.'); end
outputDir=fullfile(THIS_SCRIPT_DIR,'Data',participantID,sprintf('Visit%d',visitNumber));
if exist(outputDir,'dir')~=7, mkdir(outputDir); end
% Independent-run policy: existing files never block the experiment and are
% never overwritten. Reruns receive an EyeLink-compatible A/B/C... suffix.
fileStem=chooseAvailableRunStem(baseFileStem,outputDir);
edfFile=fileStem;
localEdfPath=fullfile(outputDir,[fileStem '.edf']);
outFile=fullfile(outputDir,[fileStem '.mat']);
emergencyFile=fullfile(outputDir,[fileStem '_INTERRUPTED.mat']);
if ~strcmp(fileStem,baseFileStem)
    warning('Existing %s output detected. This independent rerun will be saved as %s.',baseFileStem,fileStem);
end

fprintf('\n============================================================\n');
fprintf(' ER-taVNS CONTROL ASSIGNMENT\n');
fprintf('============================================================\n');
fprintf(' Participant:       %s\n',participantID);
fprintf(' Visit:             %d / 2\n',visitNumber);
fprintf(' Stimulation:       %s\n',sessionCondition);
fprintf(' Ear site:          %s\n',sessionSite);
fprintf(' Task order later:  %s -> %s\n',assignment.taskOrder{1},assignment.taskOrder{2});
fprintf(' File base name:    %s\n',fileStem);
fprintf(' Output folder:     %s\n',outputDir);
fprintf(' BrainVision name:  %s  (set manually in Recorder)\n',fileStem);
fprintf(' EyeLink/EDF name:  %s.edf\n',edfFile);
fprintf('============================================================\n\n');

%% ---------------- CONTROL PARAMETERS ------------------------------
nBursts=30;
offIntervalMinSec=P.controlOffMinSec;
offIntervalMaxSec=P.controlOffMaxSec;
rngSeed=130000+participantNumber; % same jitter sequence in both visits
rng(rngSeed,'twister');
TRIG=struct();
TRIG.STIM=10;  % same stimulation-label code used in the auditory tasks

%% ----------------------- INITIALIZE EYELINK ------------------------
% Non-fatal initialization. Missing EyeLink hardware/software is reported
% in the checklist and the control protocol continues.
eyelink_active = false;
eyelink_status_text = 'NOT CONNECTED / NOT AVAILABLE';
try
    if exist('Eyelink','file') == 0
        error('EyeLink Psychtoolbox functions are not installed/on the MATLAB path.');
    end
    if Eyelink('Initialize') ~= 0
        error('EyeLink host could not be initialized.');
    end
    Eyelink('OpenFile', edfFile, 0);
    Eyelink('command','file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    Eyelink('command','inputword_enable = YES');
    Eyelink('command','inputword_event = YES');
    Eyelink('command','inputword_sample = YES');
    Eyelink('StartRecording');
    disp('EyeLink initialized and recording.');
    eyelink_active = true;
    eyelink_status_text = 'CONNECTED';
catch ME
    errTxt = safeExceptionText(ME);
    warning('EyeLink unavailable. Continuing without EDF recording. Reason: %s', errTxt);
    eyelink_active = false;
    lowErr = lower(errTxt);
    if contains(lowErr,'invalid mex-file') || contains(lowErr,'specified module could not be found')
        eyelink_status_text = ['NOT CONNECTED / SOFTWARE DEPENDENCY ISSUE: EyeLink MEX exists but cannot load. ' ...
            'Repair/install the required EyeLink runtime/SDK or dependent Windows runtime libraries. Detail: ' errTxt];
    elseif contains(lowErr,'not installed') || contains(lowErr,'not on the matlab path')
        eyelink_status_text = ['NOT CONNECTED / SOFTWARE MISSING: ' errTxt];
    else
        eyelink_status_text = ['NOT CONNECTED / WARNING: ' errTxt];
    end
end

%% ------------------- INITIALIZE PARALLEL PORT ----------------------
% EEG / BrainVision TTL startup checklist.
%
% Simple CONNECTED criterion used here:
%   (1) Windows detects a physical parallel-port device, AND
%   (2) io64 initializes successfully and can address the configured port.
%
% This is intentionally only a startup/interface check. It does NOT insert
% any diagnostic marker into the participant recording. Actual marker values
% reaching BrainVision are verified separately beforehand with
% test_LG_trigger_codes.m / test_ERtaVNS_trigger_codes.m.
%
% On a non-lab computer without the physical LPT interface, the expected
% status is NOT CONNECTED and the experiment continues without hardware TTLs.
ioObj = [];
LPT1address = [];
eeg_status_text = 'NOT CONNECTED / NOT AVAILABLE';
eeg_ttl_output_available = false;

try
    if exist('io64','file') == 0
        error('io64 is not installed/on the MATLAB path.');
    end

    [lptState, lptDetail] = detectPhysicalParallelPort();
    % lptState: 1 = Windows detected a physical LPT device
    %           0 = Windows query succeeded and found none
    %          -1 = automatic detection unavailable/inconclusive
    if lptState == 0
        error(['Windows did not detect a physical parallel-port device. ' ...
               'The io64 software/driver alone is not treated as connected.']);
    elseif lptState < 0
        error(['Windows parallel-port detection was inconclusive. ' ...
               'EEG TTL is therefore treated as NOT CONNECTED. Detail: ' lptDetail]);
    end

    ioObj = io64;
    status = io64(ioObj);
    if status ~= 0
        error('io64 could not initialize the detected parallel-port interface.');
    end

    LPT1address = hex2dec('378');  % VERIFY in lab; PCIe cards may use another base address
    io64(ioObj, LPT1address, 0);   % put data lines in known low state

    eeg_ttl_output_available = true;
    eeg_status_text = ['CONNECTED (present physical Windows LPT detected; io64 initialized)'];

catch ME
    errTxt = safeExceptionText(ME);
    warning('EEG/parallel trigger interface unavailable. Continuing without hardware TTLs. Reason: %s', errTxt);
    ioObj = [];
    LPT1address = [];
    eeg_ttl_output_available = false;
    eeg_status_text = ['NOT CONNECTED / WARNING: ' errTxt];
end

%% ---------------------- INITIALIZE taVNS ---------------------------
stim=struct();
stim.deviceID=P.deviceID;
stim.aoChannel=P.aoChannel;
stim.gateLine=P.gateLine;
stim.hasGateLine=true;
stim.fs=P.niFs;
stim.freqHz=P.stimFrequencyHz;
stim.pulseWidthSec=P.pulseWidthSec;
stim.burstDurSec=P.stimDurationSec;
stim.preZeroSec = P.preZeroSec;
stim.postZeroSec = P.postZeroSec;
stim.taskStartLeadSec = P.taskStartLeadSec;
stim.sessionCondition=sessionCondition;
stim.sessionSite=sessionSite;
stim.syncMode='software';
[niVendorOperational,niVendorProbeText]=probeNIVendor(); %#ok<NASGU>
stim.amp=ERtaVNS_get_visit_amplitude(participantNumber,visitNumber,sessionCondition,sessionSite,true);
% The NI vendor probe is diagnostic only. It must never silently replace
% the requested stimulation amplitude with zero.
stim.niVendorProbeText=niVendorProbeText;
stim=initVNSStimulator(stim);
fprintf('TIMING CONFIG: stim %.3f s | control onset spacing %.3f-%.3f s\n',P.stimDurationSec,P.minimumStimOnsetSec,P.minimumStimOnsetSec+P.ITIJitterSec);

%% ---------------------- HARDWARE CHECKLIST -------------------------
hardwareStatus = struct();
hardwareStatus.sound = 'NOT USED BY CONTROL SCRIPT';
hardwareStatus.eyelink = eyelink_status_text;
hardwareStatus.eegBrainVision = eeg_status_text;
hardwareStatus.niDAQ = stim.niStatusText;
hardwareStatus.ds5 = stim.ds5StatusText;
hardwareStatus.sessionCondition = sessionCondition;
hardwareStatus.sessionSite = sessionSite;
hardwareStatus.vnsSyncMode = upper(stim.syncMode);
hardwareStatus.vnsAmplitude = stim.amp;
hardwareStatus.dryRun = strcmpi(stim.syncMode,'simulation');

printHardwareChecklist(hardwareStatus);

%% ---------------------- EXPERIMENT METADATA ------------------------
experimentInfo=struct();
experimentInfo.experiment='VNS_CONTROL';
experimentInfo.participantNumber=participantNumber;
experimentInfo.participantID=participantID;
experimentInfo.visitNumber=visitNumber;
experimentInfo.sessionCondition=sessionCondition;
experimentInfo.sessionSite=sessionSite;
experimentInfo.taskOrder=assignment.taskOrder;
experimentInfo.scheduleGroup=assignment.scheduleGroup;
experimentInfo.blockGroup=assignment.blockGroup;
experimentInfo.fileStem=fileStem;
experimentInfo.outputDir=outputDir;
experimentInfo.localEdfPath=localEdfPath;
experimentInfo.eyelinkHostFile=edfFile;
experimentInfo.rngSeed=rngSeed;
experimentInfo.nBursts=nBursts;
experimentInfo.offIntervalMinSec=offIntervalMinSec;
experimentInfo.offIntervalMaxSec=offIntervalMaxSec;
experimentInfo.triggerMap=TRIG;
experimentInfo.hardwareStatus=hardwareStatus;
experimentInfo.eyelink_active=eyelink_active;
stimForSave = stim;
dropStimFields = {'s','trainWave','packetAnalog','packetGate','streamData','wave','gateWave'};
presentStimFields = intersect(dropStimFields,fieldnames(stimForSave));
if ~isempty(presentStimFields)
    stimForSave = rmfield(stimForSave,presentStimFields);
end
experimentInfo.stim = stimForSave;
events_control = struct([]);
% Scalar metadata structure.  Do NOT use struct([]) here: struct([]) is a
% 0x0 structure, and MATLAB does not allow direct dot-field assignment to it.
% Initialize it as a scalar structure so the control sequence can be logged
% safely before/during the run and also in an emergency save.
controlSequenceInfo = struct();
controlSequenceInfo.condition = sessionCondition;
controlSequenceInfo.site = sessionSite;
controlSequenceInfo.offIntervalsSec = [];
sessionStartGetSecs = compatGetSecs();
experimentInfo.sessionStartGetSecs = sessionStartGetSecs;

%% ------------------------- FIXATION SCREEN -------------------------
try
    try
        fig = figure( ...
            'Color',[0.5 0.5 0.5], ...
            'MenuBar','none', ...
            'ToolBar','none', ...
            'Name','VNS control', ...
            'NumberTitle','off', ...
            'WindowState','fullscreen');
    catch
        fig = figure( ...
            'units','normalized', ...
            'outerposition',[0 0 1 1], ...
            'Color',[0.5 0.5 0.5], ...
            'MenuBar','none', ...
            'ToolBar','none', ...
            'Name','VNS control', ...
            'NumberTitle','off');
        set(fig,'Units','pixels','Position',get(0,'ScreenSize'));
    end

    axis off;
    axis equal;
    text(0.5,0.5,'Please just look at the fixation cross!', ...
        'HorizontalAlignment','center', ...
        'FontSize',48, ...
        'Color','w');
    drawnow;
    compatWaitSecs(5);

    cla;
    axis off;
    text(0.5,0.5,'+', ...
        'HorizontalAlignment','center', ...
        'FontSize',100, ...
        'Color','k');
    drawnow;
    compatWaitSecs(2);

    eventCounter=0;
    offIntervals=offIntervalMinSec + ...
        (offIntervalMaxSec-offIntervalMinSec)*rand(1,nBursts);
    controlSequenceInfo.offIntervalsSec=offIntervals;

    fprintf('\nControl: %s | %s | %d bursts\n',sessionCondition,sessionSite,nBursts);
    if eyelink_active
        Eyelink('Message',sprintf('CONTROL_BEGIN VISIT_%d CONDITION_%s SITE_%s', ...
            visitNumber,sessionCondition,sessionSite));
    end

    for burstIdx=1:nBursts
        eventCounter=eventCounter+1;
        [stim,armGetSecs]=armVNSBurst(stim);
        % Start the NI packet from a true zero baseline BEFORE the intended
        % stimulation onset.  Non-zero AO0 begins preZeroSec later.
        scheduledOnset=compatGetSecs()+1.000;
        scheduledNIStart=scheduledOnset-stim.preZeroSec;
        compatWaitSecs('UntilTime',scheduledNIStart);
        triggerGetSecs=NaN; softwareStartGetSecs=NaN;

        if strcmpi(stim.syncMode,'software')
            [stim,softwareStartGetSecs]=startVNSBurstSoftware(stim);
            triggerGetSecs=softwareStartGetSecs + stim.preZeroSec;
        elseif strcmpi(stim.syncMode,'simulation')
            softwareStartGetSecs=compatGetSecs();
            triggerGetSecs=scheduledOnset;
        end

        % Marker is placed at the intended/expected non-zero pulse onset, not
        % at the earlier NI packet start.
        compatWaitSecs('UntilTime',scheduledOnset);
        if ~isempty(ioObj)
            io64(ioObj,LPT1address,TRIG.STIM);
        end
        if eyelink_active
            Eyelink('Message',sprintf( ...
                'VNS_CONTROL BURST_%02d CODE_%d CONDITION_%s SITE_%s VISIT_%d', ...
                burstIdx,TRIG.STIM,sessionCondition,sessionSite,visitNumber));
        end
        if ~isempty(ioObj)
            compatWaitSecs(0.005);
            io64(ioObj,LPT1address,0);
        end

        thisEvent=struct( ...
            'participant',participantNumber, ...
            'visit',visitNumber, ...
            'burst',burstIdx, ...
            'condition',sessionCondition, ...
            'site',sessionSite, ...
            'marker_code',TRIG.STIM, ...
            'stim_delivery_mode',stim.syncMode, ...
            'stim_hardware_available',stim.hardwareAvailable, ...
            'arm_getsecs',armGetSecs, ...
            'scheduled_onset_getsecs',scheduledOnset, ...
            'scheduled_ni_start_getsecs',scheduledNIStart, ...
            'ni_pre_zero_ms',stim.preZeroSec*1000, ...
            'stim_trigger_getsecs',triggerGetSecs, ...
            'software_start_getsecs',softwareStartGetSecs, ...
            'trigger_schedule_offset_ms',1000*(triggerGetSecs-scheduledOnset), ...
            'burst_duration_ms',stim.burstDurSec*1000, ...
            'off_interval_after_sec',offIntervals(burstIdx));
        if isempty(events_control)
            events_control=thisEvent;
        else
            events_control(end+1)=thisEvent; %#ok<SAGROW>
        end

        nextAllowedTime=scheduledOnset+stim.burstDurSec+offIntervals(burstIdx);
        compatWaitSecs('UntilTime',nextAllowedTime);
    end

    if eyelink_active
        Eyelink('Message',sprintf('CONTROL_END VISIT_%d CONDITION_%s',visitNumber,sessionCondition));
    end

    %% ---------------------------- SAVE ------------------------------
    save(outFile,'events_control','experimentInfo','controlSequenceInfo');
    fprintf('\nControl completed.\nSaved: %s\n',outFile);
    ERtaVNS_update_tracker(participantNumber,visitNumber,'CONTROL','DONE',outFile);

catch ME
    try
        save(emergencyFile,'events_control','experimentInfo','controlSequenceInfo');
        warning('Emergency control data saved to %s', emergencyFile);
        ERtaVNS_update_tracker(participantNumber,visitNumber,'CONTROL','INTERRUPTED',emergencyFile);
    catch
    end

    try
        if ~isempty(ioObj)
            io64(ioObj,LPT1address,0);
        end
    catch
    end

    try, safelyStopVNS(stim); catch, end

    if eyelink_active
        try, Eyelink('StopRecording'); catch, end
        try, Eyelink('CloseFile'); catch, end
        try, Eyelink('ReceiveFile',edfFile,localEdfPath); catch, end
        try, Eyelink('Shutdown'); catch, end
    end

    if exist('fig','var') && ishandle(fig)
        close(fig);
    end

    rethrow(ME);
end

%% ------------------------ NORMAL SHUTDOWN --------------------------
try
    if ~isempty(ioObj)
        io64(ioObj,LPT1address,0);
    end
catch
end

try, safelyStopVNS(stim); catch, end

if eyelink_active
    try
        Eyelink('StopRecording');
        pause(1);
        Eyelink('CloseFile');
        Eyelink('ReceiveFile',edfFile,localEdfPath);
        Eyelink('Shutdown');
        disp('EyeLink control session ended and EDF retrieved.');
    catch
        warning('EyeLink shutdown/retrieval problem.');
    end
end

if exist('fig','var') && ishandle(fig)
    close(fig);
end

%% ========================== FUNCTIONS ==============================


function amp = readVNSAmplitude(label)
% Read one calibrated stimulation command value using the SAME convention
% as the working Experiment_NORA.m. No undocumented /2 conversion is applied.

    prompt = sprintf('Enter calibrated %s amplitude (Experiment_NORA convention; max 5): ', label);
    raw = input(prompt, 's');

    if contains(raw, ',')
        error('Use a dot (.) as decimal separator, not a comma (,).');
    end

    amp = str2double(raw);
    if ~isfinite(amp)
        error('A valid stimulation amplitude was not specified.');
    end
    if amp < 0 || amp > 5
        error('Stimulation amplitude must be between 0 and 5 in the lab NORA convention.');
    end
end


function stim = initVNSStimulator(stim)
% Initialize/probe the single NI AO0 -> DS5 pathway.
%
% IMPORTANT CHANGE FROM v16:
% v16 queued only the first 1-s ZERO chunk and depended on a MATLAB
% DataRequired callback to queue the later chunk containing stimulation.
% The experiment also uses blocking WaitSecs calls, so that callback can be
% delayed/starved. In that failure mode MATLAB reports IsRunning=1 while
% the NI device has only received the zero chunk.
%
% This version removes that dependency completely. For every event-related
% burst the COMPLETE short packet is queued BEFORE startBackground():
%
%   1.000 s zero -> 1.500 s 25-Hz stimulation -> 1.500 s zero
%
% The waveform generation retains the lab's Experiment_NORA convention:
%   - Dev2 / AO0
%   - port0/line0
%   - 200 kHz
%   - 25 Hz
%   - 250-us rectangular pulses
%   - amplitude entered directly in the same 0-5 command convention
%   - queueOutputData() + startBackground()
%
% The pulse cycle is derived with rectpuls exactly as in Experiment_NORA.
% It is circularly shifted so a complete pulse starts at sample 1 of the
% 1.5-s stimulation segment; therefore the intended burst onset is the
% first non-zero AO0 sample rather than ~20 ms later.

    stim.s = [];
    stim.hardwareAvailable = false;
    stim.niStatusText = 'NOT CONNECTED / NOT AVAILABLE';
    stim.ds5StatusText = 'NOT SOFTWARE-VERIFIABLE';
    stim.lastHardwareError = '';

    if ~isfield(stim,'preZeroSec'),       stim.preZeroSec = 1.000; end
    if ~isfield(stim,'postZeroSec'),      stim.postZeroSec = 1.500; end
    if ~isfield(stim,'taskStartLeadSec'), stim.taskStartLeadSec = 0.400; end

    cycleTime = 1/stim.freqHz;
    tCycle = -cycleTime/2:1/stim.fs:cycleTime/2-(1/stim.fs);
    cycleWave = stim.amp * rectpuls(tCycle,stim.pulseWidthSec);

    if stim.amp ~= 0
        firstNZ = find(cycleWave ~= 0,1,'first');
        if isempty(firstNZ)
            error('Could not construct the VNS pulse at the requested amplitude/pulse width.');
        end
        cycleWave = [cycleWave(firstNZ:end), cycleWave(1:firstNZ-1)];
    end

    nStimSamples = round(stim.burstDurSec * stim.fs);
    nCyclesNeeded = ceil(nStimSamples / numel(cycleWave));
    trainRow = repmat(cycleWave,1,nCyclesNeeded);
    trainRow = trainRow(1:nStimSamples);
    stim.trainWave = trainRow(:);

    % Force a true zero at the end of the stimulation segment.
    stim.trainWave(end) = 0;

    nPre  = round(stim.preZeroSec  * stim.fs);
    nPost = round(stim.postZeroSec * stim.fs);

    stim.packetAnalog = [zeros(nPre,1); stim.trainWave; zeros(nPost,1)];
    stim.packetGate   = [zeros(nPre,1); ones(nStimSamples,1); zeros(nPost,1)];
    stim.packetGate(end) = 0;

    stim.packetDurSec = numel(stim.packetAnalog)/stim.fs;
    stim.nPulses = nnz(diff([0; stim.trainWave~=0])==1);
    stim.actualBurstDurSec = nStimSamples/stim.fs;
    stim.actualFreqHz = stim.freqHz;
    stim.actualPulseWidthSec = stim.pulseWidthSec;
    stim.meanDutyCyclePct = 100 * stim.freqHz * stim.pulseWidthSec;

    fprintf(['VNS waveform [%s]: %.3f s train | %.3f Hz | %.1f us | ' ...
             '%d pulses | mean duty %.3f%% | pre-zero %.3f s | packet %.3f s\n'], ...
        stim.sessionCondition,stim.actualBurstDurSec,stim.actualFreqHz, ...
        stim.actualPulseWidthSec*1e6,stim.nPulses,stim.meanDutyCyclePct, ...
        stim.preZeroSec,stim.packetDurSec);
    fprintf('VNS AO0 packet QC: min %.6g V | max %.6g V | non-zero samples %d\n', ...
        min(stim.packetAnalog),max(stim.packetAnalog),nnz(stim.packetAnalog));

    % Probe the exact NI channels/API once. Do not retain this probe session.
    probeSession = [];
    try
        if exist('daqreset','file') == 0
            error('Data Acquisition Toolbox / daqreset is not available.');
        end

        daqreset;
        probeSession = daq.createSession('ni');
        probeSession.addAnalogOutputChannel(stim.deviceID,stim.aoChannel,'Voltage');

        try
            probeSession.addDigitalChannel(stim.deviceID,stim.gateLine,'OutputOnly');
            stim.hasGateLine = true;
        catch MEgate
            warning('NORA digital line unavailable; using AO0 only. Reason: %s',safeExceptionText(MEgate));
            stim.hasGateLine = false;
        end

        set(probeSession,'Rate',stim.fs);
        probeSession.IsContinuous = false;
        release(probeSession);
        probeSession = [];

        stim.hardwareAvailable = true;
        stim.niStatusText = sprintf( ...
            'CONNECTED (%s/AO%d available; COMPLETE NORA-derived packet prequeued before each burst)', ...
            stim.deviceID,stim.aoChannel);
        stim.ds5StatusText = sprintf('MANUAL CHECK REQUIRED (%s session; AO%d -> DS5)', ...
            stim.sessionCondition,stim.aoChannel);

    catch ME
        try, if ~isempty(probeSession), release(probeSession); end, catch, end
        stim.lastHardwareError = oneLineMessage(safeExceptionText(ME));
        stim.niStatusText = ['NOT CONNECTED / WARNING: ' stim.lastHardwareError];
        stim.ds5StatusText = 'NOT VERIFIABLE (NI stimulation output unavailable)';
        stim.hardwareAvailable = false;
        stim.syncMode = 'simulation';
        stim.s = [];
        warning(['NI/DS5 hardware unavailable. Continuing in SOFTWARE-SIMULATION mode. ' ...
                 'No physical stimulation will be emitted. Reason: %s'],safeExceptionText(ME));
    end

    if ~stim.hardwareAvailable
        stim.syncMode = 'simulation';
    end

    stim.isArmed = false;
    stim.lastArmGetSecs = NaN;
    stim.lastSoftwareStartCallGetSecs = NaN;
    stim.lastExpectedPulseOnsetGetSecs = NaN;
    stim.lastQueuedScans = NaN;
end


function [stim,armGetSecs] = armVNSBurst(stim)
% Prepare one event-related burst.
%
% Unlike v16, the ENTIRE packet is queued here before the NI session is
% started. Therefore stimulation delivery does not depend on DataRequired
% callbacks being serviced while Psychtoolbox WaitSecs is blocking MATLAB.

    armGetSecs = compatGetSecs();
    stim.lastArmGetSecs = armGetSecs;
    stim.isArmed = false;

    if ~isfield(stim,'hardwareAvailable') || ~stim.hardwareAvailable
        stim.syncMode = 'simulation';
        return;
    end

    try
        safelyStopVNS(stim);

        stim.s = daq.createSession('ni');
        stim.s.addAnalogOutputChannel(stim.deviceID,stim.aoChannel,'Voltage');

        if stim.hasGateLine
            stim.s.addDigitalChannel(stim.deviceID,stim.gateLine,'OutputOnly');
            outputData = [stim.packetAnalog,stim.packetGate];
        else
            outputData = stim.packetAnalog;
        end

        set(stim.s,'Rate',stim.fs);

        % Finite packet: no refill callback is required.
        stim.s.IsContinuous = false;

        % CRITICAL: queue ZERO + STIMULATION + ZERO before startBackground.
        stim.s.queueOutputData(outputData);

        stim.lastQueuedScans = size(outputData,1);
        stim.isArmed = true;

        fprintf(['VNS armed: FULL packet queued | AO0 max %.6g V | ' ...
                 'queued scans %d (%.3f s) | stimulation already present in queue\n'], ...
            max(stim.packetAnalog),stim.lastQueuedScans, ...
            stim.lastQueuedScans/stim.fs);

    catch ME
        warning('VNS arm failed; switching to simulation. Reason: %s',safeExceptionText(ME));
        stim.lastHardwareError = oneLineMessage(safeExceptionText(ME));
        stim.hardwareAvailable = false;
        stim.syncMode = 'simulation';
        try, safelyStopVNS(stim); catch, end
        stim.s = [];
    end
end


function [stim,startCallGetSecs] = startVNSBurstSoftware(stim)
% Start the already prequeued finite NORA-derived packet.
% The first non-zero AO0 sample is expected preZeroSec after this call.

    startCallGetSecs = compatGetSecs();
    stim.lastSoftwareStartCallGetSecs = startCallGetSecs;
    stim.lastExpectedPulseOnsetGetSecs = startCallGetSecs + stim.preZeroSec;

    if ~isfield(stim,'hardwareAvailable') || ~stim.hardwareAvailable || isempty(stim.s)
        stim.syncMode = 'simulation';
        return;
    end

    if ~isfield(stim,'isArmed') || ~stim.isArmed
        warning('No VNS packet was armed; treating event as simulation.');
        stim.hardwareAvailable = false;
        stim.syncMode = 'simulation';
        return;
    end

    try
        queuedBeforeStart = NaN;
        try, queuedBeforeStart = stim.s.ScansQueued; catch, end

        stim.s.startBackground();
        stim.isArmed = false;

        % This short pause is only for status reporting; the entire packet
        % was already queued, so no MATLAB callback is required for delivery.
        compatWaitSecs(0.020);

        try
            fprintf(['NI VNS packet started: IsRunning=%d | queued before start=%g | ' ...
                     'ScansQueued now=%g | expected first non-zero AO0 in %.3f s\n'], ...
                stim.s.IsRunning,queuedBeforeStart,stim.s.ScansQueued,stim.preZeroSec);
        catch
            fprintf('NI VNS packet start returned successfully; full packet was prequeued.\n');
        end

    catch ME
        warning('Software NI start failed; switching to simulation. Reason: %s',safeExceptionText(ME));
        stim.lastHardwareError = oneLineMessage(safeExceptionText(ME));
        stim.hardwareAvailable = false;
        stim.syncMode = 'simulation';
    end
end


function safelyStopVNS(stim)
% Stop/release the current finite NI packet session.

    if isempty(stim)
        return;
    end

    try
        if isfield(stim,'s') && ~isempty(stim.s)
            try
                if stim.s.IsRunning
                    stim.s.stop();
                end
            catch
            end

            try
                release(stim.s);
            catch
            end
        end
    catch
    end
end


function [ok, statusText] = probeNIVendor()
% Lightweight NI-vendor probe used before asking for stimulation amplitude.
% Never throws: it only informs dry-run behavior.
    ok = false;
    statusText = 'NI vendor not probed';
    try
        if exist('daq','file') == 0 && exist('daqreset','file') == 0
            statusText = 'SOFTWARE MISSING: Data Acquisition Toolbox functions not found';
            return;
        end

        vendors = daq.getVendors();
        if isempty(vendors)
            statusText = 'SOFTWARE/DRIVER NOT READY: no DAQ vendors found';
            return;
        end

        % Works for the table returned by current MATLAB releases.
        names = {};
        operational = [];
        try
            names = cellstr(string(vendors.ID));
        catch
            try, names = cellstr(string({vendors.ID})); catch, end
        end
        try
            operational = logical(vendors.IsOperational);
        catch
            try, operational = logical([vendors.IsOperational]); catch, end
        end

        idx = find(strcmpi(names,'ni'),1);
        if isempty(idx)
            statusText = 'SOFTWARE/DRIVER NOT READY: NI vendor support package not listed';
            return;
        end

        if isempty(operational)
            % Vendor listed but operational flag could not be parsed.
            ok = true;
            statusText = 'NI vendor listed (operational state not parsed)';
        elseif operational(idx)
            ok = true;
            statusText = 'NI vendor operational';
        else
            statusText = 'SOFTWARE/DRIVER NOT READY: NI vendor installed but not operational';
        end
    catch ME
        statusText = oneLineMessage(safeExceptionText(ME));
        ok = false;
    end
end


function printHardwareChecklist(hw)
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' HARDWARE / SOFTWARE STARTUP CHECKLIST\n');
    fprintf('============================================================\n');
    fprintf(' Sound port:              %s\n', hw.sound);
    fprintf(' EyeLink:                 %s\n', hw.eyelink);
    fprintf(' EEG / BrainVision TTL:   %s\n', hw.eegBrainVision);
    fprintf(' NI-DAQ stimulation I/O:  %s\n', hw.niDAQ);
    fprintf(' DS5 / session site:      %s\n', hw.ds5);
    fprintf(' Session condition/site:  %s / %s\n', hw.sessionCondition, hw.sessionSite);
    fprintf(' VNS synchronization:     %s\n', hw.vnsSyncMode);
    if isfield(hw,'vnsAmplitude')
        fprintf(' Requested VNS amplitude: %.6g (Experiment_NORA convention)\n', hw.vnsAmplitude);
    end

    if hw.dryRun
        fprintf('\n >>> DRY-RUN WARNING: VNS hardware is unavailable. <<<\n');
        fprintf('     Trial randomization/timing/logging will continue, but\n');
        fprintf('     NO physical taVNS/sham output will be delivered.\n');
    end

    fprintf('\nWHAT THIS CHECKLIST CAN / CANNOT TELL YOU:\n');
    fprintf(' - SOFTWARE can be checked off-lab: Psychtoolbox functions, EyeLink MEX loading,\n');
    fprintf('   io64 availability, Data Acquisition Toolbox, and NI vendor support.\n');
    fprintf(' - HARDWARE must ultimately be checked in the lab: physical EyeLink connection,\n');
    fprintf('   actual EEG TTL receipt, NI device/channel wiring, DS5 output, and final timing.\n');
    fprintf(' - EEG CONNECTED means: Windows detects a PRESENT physical LPT device and io64 initializes.\n');
    fprintf('   BrainVision marker receipt/values are checked separately with the trigger-test script.\n');
    fprintf(' - MATLAB cannot electrically identify a downstream DS5 from the NI analog output.\n');
    fprintf('============================================================\n\n');
end





function [state, detail] = detectPhysicalParallelPort()
% Conservative Windows detection of a PRESENT physical LPT interface.
%
% IMPORTANT:
%   io64 being installed/initializable is NOT enough. The experiment reports
%   EEG / BrainVision TTL = CONNECTED only when BOTH are true:
%       1) Windows currently reports a PRESENT, working physical LPT device
%          (integrated ACPI LPT or PCI/PCIe parallel card), and
%       2) io64 initializes successfully (checked by the caller).
%
% This intentionally does NOT use Win32_ParallelPort as a fallback because
% Windows can keep legacy/phantom parallel-port entries even when no usable
% physical trigger interface is present.
%
% state =  1 : present working physical LPT device detected
% state =  0 : PowerShell query succeeded; no matching present device
% state = -1 : OS query unavailable/inconclusive -> treated as NOT CONNECTED

    state = -1;
    detail = 'present physical LPT detection unavailable';

    if ~ispc
        detail = 'physical LPT auto-detection is implemented for Windows only';
        return;
    end

    try
        % Get-PnpDevice -PresentOnly is deliberately stricter than
        % Win32_ParallelPort. We accept the two common direct-I/O cases used
        % with io64: integrated ACPI parallel ports (PNP0400/PNP0401) and
        % PCI/PCIe parallel cards. A real LPT name such as (LPT1) is required.
        cmd = ['powershell -NoProfile -ExecutionPolicy Bypass -Command "' ...
               '$p=Get-PnpDevice -PresentOnly -ErrorAction SilentlyContinue | ' ...
               'Where-Object {$_.Status -eq ''OK'' -and $_.FriendlyName -match ''[(]LPT[0-9]+[)]'' -and ' ...
               '($_.InstanceId -like ''ACPI\PNP0400*'' -or $_.InstanceId -like ''ACPI\PNP0401*'' -or $_.InstanceId -like ''PCI\*'')}; ' ...
               'if($p){$p | ForEach-Object {Write-Output ($_.FriendlyName + '' | '' + $_.InstanceId + '' | '' + $_.Status)}}"'];

        [rc,out] = system(cmd);
        out = strtrim(out);

        if rc ~= 0
            detail = 'PowerShell Get-PnpDevice query failed; treating TTL interface as NOT CONNECTED';
            return;
        end

        if isempty(out)
            state = 0;
            detail = 'Windows reports no PRESENT working physical LPT device (ACPI/PCI/PCIe)';
        else
            state = 1;
            detail = ['Windows present physical LPT: ' oneLineMessage(out)];
        end
    catch ME
        detail = ['physical LPT query failed: ' oneLineMessage(safeExceptionText(ME))];
        state = -1;
    end
end


function txt = safeExceptionText(ME)
% Robust exception text WITHOUT directly accessing the exception message property.
% Some MATLAB / MEX failure paths can behave badly when code accesses
% the exception message property. getReport() is safer and preserves enough diagnostic detail.
    txt = 'unknown MATLAB error';
    try
        txt = getReport(ME, 'basic', 'hyperlinks', 'off');
    catch
        try
            txt = getReport(ME);
        catch
            txt = 'exception details unavailable';
        end
    end
    txt = oneLineMessage(txt);
end

function msg = oneLineMessage(msg)
    if isempty(msg)
        msg = 'unknown reason';
        return;
    end
    msg = regexprep(char(msg), '[\r\n\t]+', ' ');
    msg = regexprep(msg, '\s+', ' ');
    msg = strtrim(msg);
end


function t = compatGetSecs()
% Psychtoolbox GetSecs when available; otherwise a monotonic MATLAB tic/toc
% clock so the experiment logic can be tested without Psychtoolbox.
    if exist('GetSecs','file') ~= 0
        try
            t = GetSecs();
            return;
        catch
        end
    end

    persistent fallbackClock
    if isempty(fallbackClock)
        fallbackClock = tic;
    end
    t = toc(fallbackClock);
end


function compatWaitSecs(varargin)
% Psychtoolbox WaitSecs when available. MATLAB pause() fallback supports
% both WaitSecs(seconds) and WaitSecs('UntilTime', absoluteTime).
    if exist('WaitSecs','file') ~= 0
        try
            WaitSecs(varargin{:});
            return;
        catch
            % Continue to MATLAB fallback.
        end
    end

    if nargin == 1 && isnumeric(varargin{1})
        dt = max(0, double(varargin{1}));
        pause(dt);
        return;
    end

    if nargin == 2 && ischar(varargin{1}) && strcmpi(varargin{1},'UntilTime')
        target = double(varargin{2});
        dt = target - compatGetSecs();
        if dt > 0
            pause(dt);
        end
        return;
    end

    error('compatWaitSecs received unsupported arguments.');
end


function stem = chooseAvailableRunStem(baseStem,outputDir)
% Return baseStem when unused. If files already exist, append A..Z while
% respecting EyeLink's 8-character filename limit. This keeps reruns possible
% without silently overwriting prior participant/test data.
stem = baseStem;
if exist(fullfile(outputDir,[stem '.mat']),'file')~=2 && ...
        exist(fullfile(outputDir,[stem '.edf']),'file')~=2
    return;
end
letters='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
for k=1:numel(letters)
    candidate=[baseStem letters(k)];
    if length(candidate)>8
        break;
    end
    if exist(fullfile(outputDir,[candidate '.mat']),'file')~=2 && ...
            exist(fullfile(outputDir,[candidate '.edf']),'file')~=2
        stem=candidate;
        return;
    end
end
error(['No free EyeLink-compatible rerun filename is available for %s. ' ...
       'Move/archive old test files or use another participant/visit ID.'],baseStem);
end
