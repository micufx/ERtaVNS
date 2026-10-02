clc; clearvars; close all;

%% EXPERIMENT 1 (LOCAL-GLOBAL PARADIGM)
% Keep this experiment folder first on the MATLAB path (io64 files live here).
THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(THIS_SCRIPT_DIR)
    addpath(THIS_SCRIPT_DIR, '-begin');
end
try, Eyelink('Shutdown'); catch, end
clear mex
P = ERtaVNS_parameters(); % central editable stimulation/timing configuration
% v17 fix: complete NI packet is prequeued before each burst; no DataRequired refill callback.

%% ================================================================
%  EXPERIMENT 1: COUNTING | 2-SESSION EVENT-RELATED taVNS
%
%  FINAL CORE STRUCTURE (per visit):
%    - 2 sets/blocks
%    - 205 patterns per block = 25 habituation + 180 test
%    - test = 144 GS + 36 GD (80/20)
%    - per block: 18 GD stimulated + 18 GD pure-tone/no-stim
%    - 18 matched GS controls are also stimulated per block
%    - entire visit is ACTIVE (cymba) OR SHAM (earlobe)
%    - stimulation is never delivered on consecutive test trials
%    - 1850-2150 ms silence after each 650-ms series
%      -> 2500-2800 ms T5-to-T5; >=5000 ms between stimulation onsets
%
%  COUNTERBALANCING / YOKING:
%    - ACTIVE/SHAM visit order is counterbalanced across participants.
%    - E1/E2 order reverses across the two visits within participant.
%    - physical LG block order is counterbalanced across participants.
%    - LG sequence, matched-GS positions, ITI jitter, and stimulation
%      positions are yoked across ACTIVE/SHAM visits for the same task.
%
%  EEG/EYELINK MARKERS:
%    - Hardware EEG codes are deliberately limited to S1..S15.
%    - T1..T4 = S5..S8; T5 onset = S9.
%    - Local/Global condition = S1..S4 immediately after T5.
%    - Paired-stimulation label = S10 immediately after the condition marker.
%    - ACTIVE vs SHAM is fixed for the whole visit and stored in metadata.
%    - EyeLink receives the same numeric codes plus descriptive metadata.
%
%  SYNCHRONIZATION:
%    PsychPortAudio schedules sound onsets. With the currently connected
%    lab hardware, a fully prequeued NI packet is started 400 ms before T1.
%    Its 1.000-s zero pre-roll places the first non-zero AO0 pulse at T5.
%    startBackground() is therefore used without any refill callback. Before data
%    collection, quantify T5 TTL -> DS5 monitor-output latency and jitter.

%% ================================================================

%% -------------------- EXPERIMENT PARAMETERS -----------------------
fs              = P.audioFs;
toneDur         = P.toneDurSec;
soa             = P.soaSec;
patternDur      = P.patternDurSec;

% FINAL 2-SESSION TIMING:
% Silent interval from the end of the 650-ms series to the next T1 onset.
% With 1850-2150 ms silence, successive T5 onsets are 2500-2800 ms apart.
% Because paired stimulation is NEVER assigned to consecutive test trials,
% the minimum stimulation-onset-to-stimulation-onset interval is 5.0 s.
silentGapValuesSec = P.silentGapValuesSec;
seriesOnsetIntervalValuesSec = P.seriesOnsetIntervalValuesSec;

habituationTrials = 25;
nTestTrials       = 180;   % per set/block; artifact-buffered 36+36 GD design
nPatternsPerBlock = habituationTrials + nTestTrials;  % 205
nRare             = 36;    % 20% of 180
nFrequentTest     = nTestTrials - nRare;              % 144

% PILOT STIMULATION ALLOCATION PER SET/BLOCK:
%   18 GD are stimulation-paired + 18 GD remain pure-tone/no-stim.
%   18 matched GS controls are stimulation-paired.
% Across the two reversed sets this gives, per task/visit:
%   36 stimulated GD + 36 pure-tone GD + 36 stimulated matched GS.
% The +6 per GD class is an artifact buffer above the desired 30 usable
% stimulated-GD and 30 usable pure-tone-GD trials.
nStimGD = 18;
nStimGS = 18;

minStandardsBetweenRares = 2;
maxStandardsBetweenRares = 5;

amplitude       = 0.5;
rampDur         = 0.007;       % 7-ms rise/fall duration
ttlPulseDur     = 0.005;       % 5-ms TTL pulse

% Physical chords
freqs_A = [350 700 1400];
freqs_B = [500 1000 2000];

% ---------------- 4-BIT EEG TRIGGER NAMESPACE ---------------------
% Keep the numeric values identical to the BrainVision S markers expected
% on the current trigger path. The same map is used by E1 and E2.
%
%   1 = FLS condition label (Global Standard + Local Standard)
%   2 = FLD condition label (Global Standard + Local Deviant)
%   3 = RLS condition label (Global Deviant  + Local Standard)
%   4 = RLD condition label (Global Deviant  + Local Deviant)
%   5 = T1 onset
%   6 = T2 onset
%   7 = T3 onset
%   8 = T4 onset
%   9 = T5 onset (primary auditory/stimulation time anchor)
%  10 = paired-stimulation label (only on selected test trials)
%
% ACTIVE vs SHAM is not encoded by a different EEG number because the
% entire visit has one fixed stimulation condition. It is stored in the
% filename, experimentInfo, trial events, EyeLink messages and tracker.
TRIG = struct();
TRIG.FLS = 1;
TRIG.FLD = 2;
TRIG.RLS = 3;
TRIG.RLD = 4;
TRIG.T1 = 5;
TRIG.T2 = 6;
TRIG.T3 = 7;
TRIG.T4 = 8;
TRIG.T5 = 9;
TRIG.STIM_PAIRED = 10;

%% ---------------- PARTICIPANT + VISIT INFORMATION ----------------
prompt = {'Participant number / ID (e.g., S01):', 'Visit number (1 or 2):'};
dlg_title = 'ER-taVNS | Experiment 1 (Counting)';
def = {'S01','1'};
answer = inputdlg(prompt, dlg_title, 1, def);
if isempty(answer)
    error('Experiment cancelled before start.');
end

participantInput = upper(strtrim(answer{1}));
tokParticipant = regexp(participantInput, '^S?(\d+)$', 'tokens', 'once');
if isempty(tokParticipant)
    error('Participant ID must have the form S01 or a positive integer.');
end
participantNumber = str2double(tokParticipant{1});
visitNumber = str2double(strtrim(answer{2}));
if ~isfinite(participantNumber) || participantNumber < 1 || participantNumber ~= round(participantNumber)
    error('Participant number must be a positive integer.');
end
if ~ismember(visitNumber,[1 2])
    error('Visit number must be 1 or 2.');
end

[participantNumber, participantID, visitNumber] = ERtaVNS_resolve_participant(participantNumber,visitNumber,'E1',THIS_SCRIPT_DIR);
assignment = ERtaVNS_get_assignment(participantNumber,visitNumber);
sessionCondition = assignment.sessionCondition;     % ACTIVE or SHAM
blockOrder = assignment.blockOrder;                 % yoked across visits
counterbalanceGroup = assignment.blockGroup;
expectedTaskOrder = assignment.taskOrder;

% Unified EyeLink/MATLAB/BrainVision base name (<=8 EyeLink characters).
baseFileStem = sprintf('%sV%dE1',participantID,visitNumber);
if length(baseFileStem) > 8
    error('Generated base file name exceeds the EyeLink 8-character limit.');
end
outputDir = fullfile(THIS_SCRIPT_DIR,'Data',participantID,sprintf('Visit%d',visitNumber));
if exist(outputDir,'dir') ~= 7, mkdir(outputDir); end

% Scripts are intentionally runnable independently. If this component was
% already recorded, DO NOT abort and DO NOT overwrite it. Create a unique
% EyeLink-compatible rerun stem by appending A, B, C, ... (e.g. S01V1E1A).
fileStem = chooseAvailableRunStem(baseFileStem,outputDir);
edfFile = fileStem;
localEdfPath = fullfile(outputDir,[fileStem '.edf']);
outFile = fullfile(outputDir,[fileStem '.mat']);
emergencyFile = fullfile(outputDir,[fileStem '_INTERRUPTED.mat']);
if ~strcmp(fileStem,baseFileStem)
    warning('Existing %s output detected. This independent rerun will be saved as %s.',baseFileStem,fileStem);
end

fprintf('\n============================================================\n');
fprintf(' ER-taVNS ASSIGNMENT\n');
fprintf('============================================================\n');
fprintf(' Participant:       %s\n',participantID);
fprintf(' Visit:             %d / 2\n',visitNumber);
fprintf(' Stimulation:       %s\n',sessionCondition);
fprintf(' Task order today:  %s -> %s\n',expectedTaskOrder{1},expectedTaskOrder{2});
fprintf(' Running now:       E1 (Counting)\n');
fprintf(' LG block order:    %s -> %s\n',blockOrder{1},blockOrder{2});
fprintf(' File base name:    %s\n',fileStem);
fprintf(' Output folder:     %s\n',outputDir);
fprintf(' BrainVision name:  %s  (set manually in Recorder)\n',fileStem);
fprintf(' EyeLink/EDF name:  %s.edf\n',edfFile);
fprintf('============================================================\n\n');

% Warn (without aborting) when the operator starts E1 out of the assigned order.
[nextTask, trackerText] = ERtaVNS_next_task(participantNumber,visitNumber);
if ~isempty(nextTask) && ~strcmp(nextTask,'E1')
    warning('According to the tracker, the next assigned task is %s, not E1. %s',nextTask,trackerText);
end

% Tracker/order information is ADVISORY ONLY. The operator may intentionally
% run E1 at any time (for testing, restarting, or a protocol deviation).
% The warning above is retained for tracking, but it never blocks execution.

siteChoice = questdlg(sprintf(['Visit %d is %s.\nStimulation site: %s\n\n' ...
    'Confirm the correct electrode/site is connected and the visit condition is correct.'], ...
    visitNumber,sessionCondition,ERtaVNS_condition_site(sessionCondition)), ...
    'ER-taVNS site check','Confirmed','Abort','Abort');
if ~strcmp(siteChoice,'Confirmed')
    error('E1 aborted at the operator site/condition check.');
end

%% ---------------- PHYSICAL BLOCK DEFINITIONS ----------------------
% Original physical block types:
% a: AAAAA frequent / AAAAB rare
% b: BBBBB frequent / BBBBA rare
% c: AAAAB frequent / AAAAA rare
% d: BBBBA frequent / BBBBB rare
%
% The participant-specific two-block order is generated by
% ERtaVNS_get_assignment and kept IDENTICAL across Visit 1 and Visit 2.
% This yokes physical context, random LG sequence, and stimulation positions
% between ACTIVE and SHAM visits.
blockDefs = struct();
blockDefs.a.name = 'a'; blockDefs.a.frequentPattern = 'AAAAA'; blockDefs.a.rarePattern = 'AAAAB';
blockDefs.b.name = 'b'; blockDefs.b.frequentPattern = 'BBBBB'; blockDefs.b.rarePattern = 'BBBBA';
blockDefs.c.name = 'c'; blockDefs.c.frequentPattern = 'AAAAB'; blockDefs.c.rarePattern = 'AAAAA';
blockDefs.d.name = 'd'; blockDefs.d.frequentPattern = 'BBBBA'; blockDefs.d.rarePattern = 'BBBBB';

%% ---------------- REPRODUCIBLE / YOKED RANDOMIZATION -------------
% IMPORTANT: visit number and stimulation condition are deliberately NOT in
% this seed. Thus the same participant receives the same LG trial order,
% matched-GS positions, ITI jitter, and stimulated-trial positions in both
% visits for E1. This minimizes period/session confounding.
rngSeed = 110000 + participantNumber;
rng(rngSeed,'twister');

%% -------------------- BUILD PHYSICAL TONES ------------------------
toneA = makeChord(freqs_A, fs, toneDur, amplitude, rampDur);
toneB = makeChord(freqs_B, fs, toneDur, amplitude, rampDur);

%% --------------------- INITIALIZE AUDIO ----------------------------
% Non-fatal initialization: the experiment can be dry-run on a computer
% without an audio device / PsychPortAudio. When audio is unavailable, all
% trial timing, randomization, figures and event logs still run, but no
% acoustic signal is physically played.
audio_active = false;
pahandle = [];
audio_status_text = 'NOT CONNECTED / NOT AVAILABLE';

try
    if exist('InitializePsychSound','file') == 0 || exist('PsychPortAudio','file') == 0
        error('Psychtoolbox PsychPortAudio functions are not installed/on the MATLAB path.');
    end
    InitializePsychSound(1);

    % Two output channels; the same signal is sent to left and right.
    % Latency class 1 is a conservative low-latency setting.
    pahandle = PsychPortAudio('Open', [], 1, 1, fs, 2);
    PsychPortAudio('Volume', pahandle, 1);
    audio_active = true;
    audio_status_text = 'CONNECTED';
catch ME
    warning('Sound output unavailable. Continuing in silent software-test mode. Reason: %s', safeExceptionText(ME));
    audio_active = false;
    pahandle = [];
    audio_status_text = ['NOT CONNECTED / WARNING: ' oneLineMessage(safeExceptionText(ME))];
end

%% ----------------------- INITIALIZE EYELINK ------------------------
% Non-fatal initialization. Missing EyeLink hardware/software is reported
% in the checklist and the rest of the experiment continues.
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
% FINAL 2-SESSION HARDWARE ARCHITECTURE:
% One DS5 is used. The physical ear site is changed BETWEEN visits:
%   ACTIVE visit -> left cymba conchae
%   SHAM visit   -> left earlobe
% Therefore only the already-used NI Dev2/AO0 channel is required.
%
% With the hardware currently connected in the lab, NI stimulation starts
% with a fully prequeued NORA-derived NI packet started by startBackground(). The exact
% T5-to-current latency/jitter MUST be bench-measured from the DS5 monitor
% output before participant acquisition. A hardware-PFI mode can be added
% later if a dedicated trigger cable is installed.
stim = struct();
stim.deviceID = P.deviceID;
stim.aoChannel = P.aoChannel;                    % existing AO0 -> DS5 rear BNC INPUT
stim.gateLine = P.gateLine;         % preserved from Experiment_NORA
stim.hasGateLine = true;

stim.fs = P.niFs;
stim.freqHz = P.stimFrequencyHz;
stim.pulseWidthSec = P.pulseWidthSec;
stim.burstDurSec = P.stimDurationSec;
stim.preZeroSec = P.preZeroSec; % derived so non-zero onset remains at T5
stim.postZeroSec = P.postZeroSec;
stim.taskStartLeadSec = P.taskStartLeadSec;
stim.sessionCondition = sessionCondition;
stim.sessionSite = ERtaVNS_condition_site(sessionCondition);
stim.syncMode = 'software';            % current connected lab configuration

[niVendorOperational, niVendorProbeText] = probeNIVendor(); %#ok<NASGU>
stim.amp = ERtaVNS_get_visit_amplitude(participantNumber,visitNumber,sessionCondition,stim.sessionSite,false);
% The NI vendor probe is diagnostic only. It must never silently replace
% the requested stimulation amplitude with zero.
stim.niVendorProbeText = niVendorProbeText;
stim = initVNSStimulator(stim);
fprintf(['TIMING CONFIG: stim %.3f s | %.1f Hz | %.0f us | silent gap %.0f-%.0f ms | ' ...
         'earliest stim-onset spacing %.3f s\n'], ...
    P.stimDurationSec,P.stimFrequencyHz,P.pulseWidthSec*1e6, ...
    1000*P.silentGapMinSec,1000*P.silentGapMaxSec,P.actualMinimumStimOnsetSec);

%% ---------------------- HARDWARE CHECKLIST -------------------------
hardwareStatus = struct();
hardwareStatus.sound = audio_status_text;
hardwareStatus.eyelink = eyelink_status_text;
hardwareStatus.eegBrainVision = eeg_status_text;
hardwareStatus.niDAQ = stim.niStatusText;
hardwareStatus.ds5 = stim.ds5StatusText;
hardwareStatus.sessionCondition = sessionCondition;
hardwareStatus.sessionSite = stim.sessionSite;
hardwareStatus.vnsSyncMode = upper(stim.syncMode);
hardwareStatus.vnsAmplitude = stim.amp;
hardwareStatus.dryRun = strcmpi(stim.syncMode,'simulation');

printHardwareChecklist(hardwareStatus);

%% ---------------------- EXPERIMENT METADATA ------------------------
experimentInfo = struct();
experimentInfo.experiment = 1;
experimentInfo.participantNumber = participantNumber;
experimentInfo.participantID = participantID;
experimentInfo.fileStem = fileStem;
experimentInfo.outputDir = outputDir;
experimentInfo.localEdfPath = localEdfPath;
experimentInfo.eyelinkHostFile = edfFile;
experimentInfo.visitNumber = visitNumber;
experimentInfo.sessionCondition = sessionCondition;
experimentInfo.sessionSite = stim.sessionSite;
experimentInfo.taskOrder = expectedTaskOrder;
experimentInfo.scheduleGroup = assignment.scheduleGroup;
experimentInfo.counterbalanceGroup = counterbalanceGroup;
experimentInfo.blockOrder = blockOrder;
experimentInfo.rngSeed = rngSeed;

experimentInfo.fs = fs;
experimentInfo.toneDur = toneDur;
experimentInfo.soa = soa;
experimentInfo.patternDur = patternDur;
experimentInfo.silentGapValuesSec = silentGapValuesSec;
experimentInfo.silentGapMinSec = min(silentGapValuesSec);
experimentInfo.silentGapMaxSec = max(silentGapValuesSec);
experimentInfo.silentGapStepSec = P.ITIStepSec;
experimentInfo.silentGapDefinition = 'silence from end of 650-ms series to onset of next series';
experimentInfo.seriesOnsetIntervalValuesSec = seriesOnsetIntervalValuesSec;
experimentInfo.seriesOnsetIntervalMinSec = min(seriesOnsetIntervalValuesSec);
experimentInfo.seriesOnsetIntervalMaxSec = max(seriesOnsetIntervalValuesSec);
experimentInfo.rampDur = rampDur;
experimentInfo.ttlPulseDur = ttlPulseDur;

experimentInfo.nStimGD_perBlock = nStimGD;
experimentInfo.nPureToneGD_perBlock = nRare - nStimGD;
experimentInfo.nStimGS_perBlock = nStimGS;
experimentInfo.targetUsableStimGD_perTaskVisit = 30;
experimentInfo.targetUsablePureToneGD_perTaskVisit = 30;
experimentInfo.artifactBufferPerGDClass = 6;
experimentInfo.noConsecutiveStimulation = true;
experimentInfo.minimumStimOnsetIntervalSec = P.minimumStimOnsetSec;
experimentInfo.actualMinimumStimOnsetIntervalSec = 2 * min(seriesOnsetIntervalValuesSec);

stimForSave = stim;
dropStimFields = {'s','trainWave','packetAnalog','packetGate','streamData','wave','gateWave'};
presentStimFields = intersect(dropStimFields,fieldnames(stimForSave));
if ~isempty(presentStimFields)
    stimForSave = rmfield(stimForSave,presentStimFields);
end
experimentInfo.stim = stimForSave;
experimentInfo.freqs_A = freqs_A;
experimentInfo.freqs_B = freqs_B;

experimentInfo.habituationTrials = habituationTrials;
experimentInfo.nPatternsPerBlock = nPatternsPerBlock;
experimentInfo.nTestTrials = nTestTrials;
experimentInfo.nRare = nRare;
experimentInfo.nFrequentTest = nFrequentTest;
experimentInfo.minStandardsBetweenRares = minStandardsBetweenRares;
experimentInfo.maxStandardsBetweenRares = maxStandardsBetweenRares;
experimentInfo.triggerMap = TRIG;
experimentInfo.hardwareStatus = hardwareStatus;
experimentInfo.audio_active = audio_active;
experimentInfo.eyelink_active = eyelink_active;

events_audio_01 = struct([]);
blockSequenceInfo = struct([]);
urevent_counter = 1;

% Command-window monitoring only. This counter has no role in timing,
% triggering, stimulation delivery, or saved trial classification.
stimConsoleCounter = 0;

% One continuous clock for the entire two-block session
sessionStartGetSecs = compatGetSecs();
experimentInfo.sessionStartGetSecs = sessionStartGetSecs;

%% =========================== RUN BLOCKS ============================
try
    for presentedBlock = 1:2

        blockKey = blockOrder{presentedBlock};
        binfo = blockDefs.(blockKey);

        fprintf('\n============================================\n');
        fprintf('Presented block %d/2: type %s\n', presentedBlock, blockKey);
        fprintf('Frequent: %s | Rare: %s\n', ...
            binfo.frequentPattern, binfo.rarePattern);
        fprintf('============================================\n');

        %% ---------- Construct constrained 180-pattern test phase ----
        [testIsRare, standardsBeforeEachRare, trailingStandards] = ...
            makeConstrainedRareSequence( ...
                nFrequentTest, nRare, ...
                minStandardsBetweenRares, ...
                maxStandardsBetweenRares);

        if length(testIsRare) ~= nTestTrials
            error('Internal sequence-generation error: test phase length is incorrect.');
        end
        if sum(testIsRare) ~= nRare
            error('Internal sequence-generation error: rare count is incorrect.');
        end

        % 25 habituation patterns + 180 test patterns
        patternList = cell(1, nPatternsPerBlock);
        phaseList   = cell(1, nPatternsPerBlock);
        isRareList  = false(1, nPatternsPerBlock);

        for k = 1:habituationTrials
            patternList{k} = binfo.frequentPattern;
            phaseList{k}   = 'Habituation';
            isRareList(k)  = false;
        end

        for k = 1:nTestTrials
            idx = habituationTrials + k;
            phaseList{idx} = 'Test';
            isRareList(idx) = testIsRare(k);

            if testIsRare(k)
                patternList{idx} = binfo.rarePattern;
            else
                patternList{idx} = binfo.frequentPattern;
            end
        end

        % Save exact sequence information
        blockSequenceInfo(presentedBlock).presentedBlock = presentedBlock;
        blockSequenceInfo(presentedBlock).blockType = blockKey;
        blockSequenceInfo(presentedBlock).frequentPattern = binfo.frequentPattern;
        blockSequenceInfo(presentedBlock).rarePattern = binfo.rarePattern;
        blockSequenceInfo(presentedBlock).testIsRare = testIsRare;
        blockSequenceInfo(presentedBlock).standardsBeforeEachRare = standardsBeforeEachRare;
        blockSequenceInfo(presentedBlock).trailingStandards = trailingStandards;
        blockSequenceInfo(presentedBlock).fullPatternList = patternList;

        %% ---------- Balanced event-related VNS assignment ------------
        % Habituation is NEVER stimulated.
        % In each 180-trial test phase:
        %   - 18 of 36 GD are stimulated with the VISIT condition;
        %   - the other 18 GD are pure-tone / unstimulated;
        %   - 18 matched GS controls are stimulated with the same condition;
        %   - the remaining 126 GS are unstimulated;
        %   - stimulation is NEVER assigned on consecutive test trials.
        % The exact positions are deterministic/yoked across Visit 1 and 2.
        [testStimPlan, matchedGSIdxTest, stimPlanQC] = makeSessionVNSPlan( ...
            testIsRare, nStimGD, nStimGS, sessionCondition);

        fullStimPlan = [repmat({'NONE'},1,habituationTrials), testStimPlan];
        fullPostStimBuffer = [false(1,habituationTrials), stimPlanQC.postStimBufferMask];

        blockSequenceInfo(presentedBlock).testStimPlan = testStimPlan;
        blockSequenceInfo(presentedBlock).fullStimPlan = fullStimPlan;
        blockSequenceInfo(presentedBlock).stimulatedGDIdxTest = stimPlanQC.stimulatedGDIdx;
        blockSequenceInfo(presentedBlock).pureToneGDIdxTest = stimPlanQC.pureToneGDIdx;
        blockSequenceInfo(presentedBlock).matchedGSIdxTest = matchedGSIdxTest;
        blockSequenceInfo(presentedBlock).stimPlanQC = stimPlanQC;
        blockSequenceInfo(presentedBlock).postStimBufferIdxTest = stimPlanQC.postStimBufferIdx;
        blockSequenceInfo(presentedBlock).fullPostStimBuffer = fullPostStimBuffer;

        fprintf(['Block %d VNS plan [%s]: GD=%d, matched GS=%d, stimulated total=%d, ' ...
                 'post-stim GS buffers=%d, pure-tone GD=%d, minimum spacing=%d trial(s)\n'], ...
            presentedBlock,sessionCondition,stimPlanQC.nGDStimulated,stimPlanQC.nGS, ...
            stimPlanQC.nStimulated,stimPlanQC.nPostStimBuffers,stimPlanQC.nGDPureTone, ...
            stimPlanQC.minTrialSeparation);

        %% ---------- Silent inter-series interval ----------------------
        % Final protocol: 1850-2150 ms silence after each 650-ms series.
        % With no consecutive stimulation this guarantees >=5.0 s between
        % possible stimulation onsets.
        silentGapSec = silentGapValuesSec(randi(numel(silentGapValuesSec), ...
            1, nPatternsPerBlock-1));
        seriesOnsetIntervalSec = patternDur + silentGapSec;

        blockSequenceInfo(presentedBlock).silentGapSec = silentGapSec;
        blockSequenceInfo(presentedBlock).silentGap_ms = silentGapSec * 1000;
        blockSequenceInfo(presentedBlock).seriesOnsetIntervalSec = seriesOnsetIntervalSec;
        blockSequenceInfo(presentedBlock).seriesOnsetInterval_ms = seriesOnsetIntervalSec * 1000;

        %% -------------------- Fixation window ------------------------
        try
            fig = figure( ...
                'Color',[0.5 0.5 0.5], ...
                'MenuBar','none', ...
                'ToolBar','none', ...
                'Name',sprintf('Block %d', presentedBlock), ...
                'NumberTitle','off', ...
                'WindowState','fullscreen');
        catch
            fig = figure( ...
                'units','normalized', ...
                'outerposition',[0 0 1 1], ...
                'Color',[0.5 0.5 0.5], ...
                'MenuBar','none', ...
                'ToolBar','none', ...
                'Name',sprintf('Block %d', presentedBlock), ...
                'NumberTitle','off');
            set(fig,'Units','pixels','Position',get(0,'ScreenSize'));
        end

        axis off;
        axis equal;

        if presentedBlock == 1
            text(0.5, 0.58, 'Please focus on the sounds', ...
                'HorizontalAlignment','center', ...
                'FontSize',48, ...
                'Color','w');
            text(0.5, 0.42, 'and silently count all rare (deviant) sound patterns.', ...
                'HorizontalAlignment','center', ...
                'FontSize',48, ...
                'Color','w');
            drawnow;
            compatWaitSecs(5);
            cla;
            axis off;
        end

        text(0.5, 0.5, '+', ...
            'HorizontalAlignment','center', ...
            'FontSize',100, ...
            'Color','k');
        drawnow;

        % Give participant 2 s of fixation before first pattern.
        compatWaitSecs(2);

        % Schedule first pattern slightly in the future. Subsequent pattern
        % onsets follow the pre-generated 1850-2150 ms SILENT-gap sequence.
        nextPatternOnset = compatGetSecs() + 0.500;

        %% ======================= PATTERN LOOP ========================
        for t = 1:nPatternsPerBlock

            pattern = patternList{t};
            phaseName = phaseList{t};
            isHabituation = strcmp(phaseName, 'Habituation');

            % Timing from THIS series to the next one.
            % silent gap = end of current 650-ms series -> next T1 onset
            % onset interval = current T1 onset -> next T1 onset
            if t < nPatternsPerBlock
                thisSilentGapSec = silentGapSec(t);
                thisSilentGap_ms = thisSilentGapSec * 1000;
                thisSeriesOnsetIntervalSec = seriesOnsetIntervalSec(t);
                thisSeriesOnsetInterval_ms = thisSeriesOnsetIntervalSec * 1000;
            else
                thisSilentGapSec = NaN;
                thisSilentGap_ms = NaN;
                thisSeriesOnsetIntervalSec = NaN;
                thisSeriesOnsetInterval_ms = NaN;
            end

            stimCondition = fullStimPlan{t};
            isStimulatedTrial = ~strcmp(stimCondition,'NONE');
            isPostStimBuffer = fullPostStimBuffer(t);

            % ---------------- Determine local/global labels ----------
            if pattern(1) == pattern(5)
                localLabel = 'Local Standard';
            else
                localLabel = 'Local Deviant';
            end

            if strcmp(pattern, binfo.frequentPattern)
                globalLabel = 'Global Standard';
            else
                globalLabel = 'Global Deviant';
            end

            % Full condition is known BEFORE the pattern is played because
            % both the block context and the trial sequence are predetermined.
            % Therefore the fifth-tone TTL can carry the correct condition
            % exactly at fifth-tone onset; no post-hoc block inference is needed.
            if strcmp(localLabel,'Local Standard') && strcmp(globalLabel,'Global Standard')
                typeLabel = 'Frequent Local Standard';
                typeShort = 'FLS';
                localShort = 'LS';
                globalShort = 'GS';
                conditionIndex = 1;
            elseif strcmp(localLabel,'Local Deviant') && strcmp(globalLabel,'Global Standard')
                typeLabel = 'Frequent Local Deviant';
                typeShort = 'FLD';
                localShort = 'LD';
                globalShort = 'GS';
                conditionIndex = 2;
            elseif strcmp(localLabel,'Local Standard') && strcmp(globalLabel,'Global Deviant')
                typeLabel = 'Rare Local Standard';
                typeShort = 'RLS';
                localShort = 'LS';
                globalShort = 'GD';
                conditionIndex = 3;
            elseif strcmp(localLabel,'Local Deviant') && strcmp(globalLabel,'Global Deviant')
                typeLabel = 'Rare Local Deviant';
                typeShort = 'RLD';
                localShort = 'LD';
                globalShort = 'GD';
                conditionIndex = 4;
            else
                error('Unexpected local/global condition.');
            end
            % ---------------- T5 + condition/stimulation markers ---------
            % S9 is written at the scheduled fifth-tone onset and is the
            % primary EEG anchor for both the fifth tone and VNS onset.
            % After S9 returns low, S1..S4 identify the LG condition and,
            % on selected test trials, S10 labels that stimulation was paired.
            fifthPhysical = pattern(5);
            t5Code = TRIG.T5;
            t5Label = typeShort;
            trialConditionCode = TRIG.(typeShort);

            sessionTrial = (presentedBlock-1)*nPatternsPerBlock + t;

            % ---------------- Arm VNS burst before the pattern --------
            stimArmGetSecs = NaN;
            stimTriggerGetSecs = NaN;
            stimSoftwareStartCallGetSecs = NaN;
            stimMarkerCode = NaN;
            if isStimulatedTrial
                [stim, stimArmGetSecs] = armVNSBurst(stim);
                stimMarkerCode = TRIG.STIM_PAIRED;
            end

            % ---------------- Make complete 650-ms waveform ----------
            monoWave = makePatternWaveform(pattern, toneA, toneB, fs, soa, toneDur);
            stereoWave = [monoWave; monoWave];
            scheduledPatternOnset = nextPatternOnset;

            if audio_active
                PsychPortAudio('FillBuffer', pahandle, stereoWave);
                PsychPortAudio('Start', pahandle, 1, scheduledPatternOnset, 0);
            end

            % NORA-derived finite packet: the NI packet begins 400 ms before
            % T1 with 1000 ms of zero AO0, so its first non-zero pulse is
            % expected at T5 (= T1 + 600 ms).
            scheduledNIStartGetSecs = NaN;
            if isStimulatedTrial
                scheduledNIStartGetSecs = scheduledPatternOnset - stim.taskStartLeadSec;
                compatWaitSecs('UntilTime',scheduledNIStartGetSecs);
                if strcmpi(stim.syncMode,'software')
                    [stim, stimSoftwareStartCallGetSecs] = startVNSBurstSoftware(stim);
                    stimTriggerGetSecs = stimSoftwareStartCallGetSecs + stim.preZeroSec;
                elseif strcmpi(stim.syncMode,'simulation')
                    stimSoftwareStartCallGetSecs = compatGetSecs();
                    stimTriggerGetSecs = scheduledPatternOnset + 4*soa;
                end
            end

            % ------------------ Per-tone onset triggers ---------------
            % Hardware/BrainVision codes are deliberately <=15:
            % T1..T4 = 5,6,7,8 and T5 = 9. After T5, S1..S4 label
            % LG condition and S10 labels paired stimulation. EyeLink mirrors
            % the same numeric codes in descriptive messages.
            toneTriggerCodes = nan(1,5);
            ttlWriteGetSecs = NaN;

            for s = 1:5
                scheduledToneOnset = scheduledPatternOnset + (s-1)*soa;
                compatWaitSecs('UntilTime', scheduledToneOnset);
                physicalSound = pattern(s);

                if s == 1
                    toneCode = TRIG.T1;
                elseif s == 2
                    toneCode = TRIG.T2;
                elseif s == 3
                    toneCode = TRIG.T3;
                elseif s == 4
                    toneCode = TRIG.T4;
                else
                    toneCode = t5Code;
                end
                toneTriggerCodes(s) = toneCode;

                if ~isempty(ioObj)
                    io64(ioObj, LPT1address, toneCode);
                    ttlWriteGetSecs = compatGetSecs();
                end
                % NI stream was already started before T1. Its 1-s zero pre-roll
                % places the expected first non-zero AO0 pulse at T5.

                if eyelink_active
                    if s <= 4
                        Eyelink('Message',sprintf( ...
                            'LG_TONE TRIAL_%03d BLOCK_%d POS_%d CODE_%d PHYS_%c', ...
                            sessionTrial,presentedBlock,s,toneCode,physicalSound));
                    else
                        Eyelink('Message',sprintf( ...
                            ['LG_T5 TRIAL_%03d BLOCK_%d POS_5 CODE_%d PHASE_%s ' ...
                             'TYPE_%s LOCAL_%s GLOBAL_%s PHYS_%c PATTERN_%s ' ...
                             'STIM_%d SESSION_%s VISIT_%d'], ...
                            sessionTrial,presentedBlock,toneCode,upper(phaseName), ...
                            t5Label,localShort,globalShort,physicalSound,pattern, ...
                            isStimulatedTrial,sessionCondition,visitNumber));
                    end
                end

                % Hold each auditory marker for 5 ms, then return data lines low.
                if ~isempty(ioObj)
                    compatWaitSecs(ttlPulseDur);
                    io64(ioObj,LPT1address,0);
                end

                % After the T5 onset marker, emit descriptive LG/stimulation
                % labels. Physical stimulation has already been started at T5.
                if s == 5
                    compatWaitSecs(0.001);
                    if ~isempty(ioObj)
                        io64(ioObj,LPT1address,trialConditionCode);
                        compatWaitSecs(ttlPulseDur);
                        io64(ioObj,LPT1address,0);
                    end
                    if eyelink_active
                        Eyelink('Message',sprintf( ...
                            'LG_CONDITION TRIAL_%03d BLOCK_%d CODE_%d TYPE_%s', ...
                            sessionTrial,presentedBlock,trialConditionCode,typeShort));
                    end
                    % Buffer identity is NOT sent as a dedicated EEG TTL or EyeLink message.
                    % It is reconstructed during preprocessing as trial N+1 after a
                    % stimulation-paired trial and verified to be an unstimulated GS.

                    if isStimulatedTrial
                        compatWaitSecs(0.001);
                        if ~isempty(ioObj)
                            io64(ioObj,LPT1address,stimMarkerCode);
                            compatWaitSecs(ttlPulseDur);
                            io64(ioObj,LPT1address,0);
                        end
                        if eyelink_active
                            Eyelink('Message',sprintf( ...
                                'VNS_TRIAL TRIAL_%03d BLOCK_%d CODE_%d CONDITION_%s SITE_%s', ...
                                sessionTrial,presentedBlock,stimMarkerCode,sessionCondition,stim.sessionSite));
                        end

                        % Non-time-critical command-window monitoring. Print only
                        % AFTER the T5, condition and S10 marker operations above, so
                        % console I/O cannot delay the scheduled T5/VNS onset.
                        stimConsoleCounter = stimConsoleCounter + 1;
                        fprintf('[VNS] #%02d | Block %d | Trial %03d | %s | %s | next trial = post-stim GS buffer\n', ...
                            stimConsoleCounter, presentedBlock, sessionTrial, globalShort, sessionCondition);
                    end
                end

            end

            % Wait until physical pattern has ended
            scheduledPatternEnd = scheduledPatternOnset + patternDur;
            compatWaitSecs('UntilTime', scheduledPatternEnd);

            % Retrieve PsychPortAudio's reported actual device start time
            % only when a real audio device is active.
            actualPatternOnset = NaN;
            if audio_active
                try
                    audioStatus = PsychPortAudio('GetStatus', pahandle);
                    if isfield(audioStatus,'StartTime') && audioStatus.StartTime > 0
                        actualPatternOnset = audioStatus.StartTime;
                    end
                catch ME
                    warning('PsychPortAudio status failed; continuing with scheduled timing. Reason: %s', safeExceptionText(ME));
                end

                % Ensure playback is finished/stopped before refilling next trial.
                try
                    PsychPortAudio('Stop', pahandle, 1);
                catch
                end
            end

            % LG condition/stimulation labels were emitted immediately after the T5
            % onset marker; there is no additional post-pattern trigger.

            % ---------------------- Timing log -----------------------
            scheduledFifthToneOnset = scheduledPatternOnset + 4*soa;

            if isfinite(actualPatternOnset)
                actualFifthToneOnset = actualPatternOnset + 4*soa;
                timingReference = actualFifthToneOnset;
            else
                actualFifthToneOnset = NaN;
                timingReference = scheduledFifthToneOnset;
            end

            sessionTimeSec = timingReference - sessionStartGetSecs;
            sessionTimeMs  = sessionTimeSec * 1000;

            % Correct conversion from seconds to samples.
            % This is an AUDIO-clock sample index, not an EEG sample index.
            latencyAudioSamples = round(sessionTimeSec * fs) + 1;

            if strcmp(localShort,'LS')
                structuralPattern = [1 1 1 1 1];
            else
                structuralPattern = [1 1 1 1 2];
            end

            if isHabituation
                modalityLabel = 'habituation';
            else
                modalityLabel = 'test';
            end

            newEvent = struct( ...
                'type',                         typeLabel, ...
                'effect',                       globalLabel, ...
                'sound',                        localLabel, ...
                'experiment',                   1, ...
                'pattern',                      structuralPattern, ...
                'modality',                     modalityLabel, ...
                'block',                        presentedBlock, ...
                'latency',                      latencyAudioSamples, ...
                'urevent',                      urevent_counter, ...
                'duration',                     1, ...
                'time',                         sessionTimeMs, ...
                'physical_pattern',             pattern, ...
                'phase',                        phaseName, ...
                'is_habituation',               isHabituation, ...
                'is_rare_test',                 (~isHabituation && isRareList(t)), ...
                'participant',                  participantNumber, ...
                'visit',                        visitNumber, ...
                'session_condition',            sessionCondition, ...
                'counterbalance_group',         counterbalanceGroup, ...
                'presented_block',              presentedBlock, ...
                'physical_block_type',          blockKey, ...
                'trial_in_block',               t, ...
                'session_trial',                sessionTrial, ...
                'silent_gap_to_next_ms',         thisSilentGap_ms, ...
                'series_onset_to_next_ms',      thisSeriesOnsetInterval_ms, ...
                'series_ITI_to_next_ms',        thisSeriesOnsetInterval_ms, ... % legacy field name
                'type_short',                   typeShort, ...
                'local_short',                  localShort, ...
                'global_short',                 globalShort, ...
                'condition_index',              conditionIndex, ...
                'is_post_stim_buffer',          isPostStimBuffer, ...
                'exclude_primary_eeg',          isPostStimBuffer, ...
                'exclude_primary_pupil',        isPostStimBuffer, ...
                'stim_condition',               stimCondition, ...
                'stim_site',                    stimConditionToSite(stimCondition), ...
                'stim_delivery_mode',           stim.syncMode, ...
                'stim_hardware_available',      stim.hardwareAvailable, ...
                'is_stimulated',                isStimulatedTrial, ...
                'is_matched_gs_control',        (~isHabituation && strcmp(globalShort,'GS') && isStimulatedTrial), ...
                'is_stimulated_gd',             (~isHabituation && strcmp(globalShort,'GD') && isStimulatedTrial), ...
                'is_pure_tone_gd',              (~isHabituation && strcmp(globalShort,'GD') && ~isStimulatedTrial), ...
                'stim_marker_code',             stimMarkerCode, ...
                'stim_arm_getsecs',             stimArmGetSecs, ...
                'stim_trigger_getsecs',         stimTriggerGetSecs, ...
                'stim_software_start_getsecs',  stimSoftwareStartCallGetSecs, ...
                'stim_expected_burst_ms',       stim.burstDurSec*1000, ...
                'stim_t5_offset_ms',            1000*(stimTriggerGetSecs - scheduledFifthToneOnset), ...
                'stim_actual_t5_offset_ms',     1000*(stimTriggerGetSecs - actualFifthToneOnset), ...
                't5_ttl_write_getsecs',         ttlWriteGetSecs, ...
                't5_trigger_code',              t5Code, ...
                't5_trigger_label',             t5Label, ...
                'trial_condition_code',         trialConditionCode, ...
                'tone_trigger_codes',           toneTriggerCodes, ...
                'fifth_tone_physical',          pattern(5), ...
                'scheduled_pattern_onset',      scheduledPatternOnset, ...
                'actual_pattern_onset',         actualPatternOnset, ...
                'scheduled_fifth_tone_onset',   scheduledFifthToneOnset, ...
                'actual_fifth_tone_onset',      actualFifthToneOnset, ...
                'time_sec',                     sessionTimeSec, ...
                'time_ms',                      sessionTimeMs, ...
                'latency_audio_samples',        latencyAudioSamples);

            if isempty(events_audio_01)
                events_audio_01 = newEvent;
            else
                events_audio_01(end+1) = newEvent; %#ok<SAGROW>
            end

            urevent_counter = urevent_counter + 1;

            % Schedule the next pattern using the Bekinschtein-style
            % 1850-2150 ms silent gap after the series (50-ms steps).
            if t < nPatternsPerBlock
                nextPatternOnset = scheduledPatternOnset + thisSeriesOnsetIntervalSec;
            end
        end

        if ishandle(fig)
            close(fig);
        end

        %% ----------------------- Between blocks ---------------------
        if presentedBlock == 1
            choice = questdlg( ...
                sprintf(['Block 1 is complete.\n\n' ...
                         'Next block type: %s\n\nContinue?'], blockOrder{2}), ...
                'Next Block', ...
                'Yes','No','Yes');

            if strcmp(choice,'No')
                warning('Experiment stopped by operator after Block 1.');
                break;
            end
        end
    end

    %% -------------------------- SAVE DATA ---------------------------
    save(outFile, ...
        'events_audio_01', ...
        'experimentInfo', ...
        'blockSequenceInfo');

    fprintf('\nExperiment completed.\nSaved: %s\n', outFile);
    ERtaVNS_update_tracker(participantNumber,visitNumber,'E1','DONE',outFile);

catch ME
    % Save whatever was acquired before stopping.
    try
        save(emergencyFile, ...
            'events_audio_01', ...
            'experimentInfo', ...
            'blockSequenceInfo');
        warning('Emergency data saved to %s', emergencyFile);
        ERtaVNS_update_tracker(participantNumber,visitNumber,'E1','INTERRUPTED',emergencyFile);
    catch
    end

    try
        if ~isempty(ioObj)
            io64(ioObj, LPT1address, 0);
        end
    catch
    end

    try
        if exist('stim','var')
            safelyStopVNS(stim);
        end
    catch
    end

    if audio_active && ~isempty(pahandle), try, PsychPortAudio('Close', pahandle); catch, end, end

    if eyelink_active
        try, Eyelink('StopRecording'); catch, end
        try, Eyelink('CloseFile'); catch, end
        try, Eyelink('ReceiveFile', edfFile, localEdfPath); catch, end
        try, Eyelink('Shutdown'); catch, end
    end

    if exist('fig','var') && ishandle(fig)
        close(fig);
    end

    fprintf(2,'\nORIGINAL ERROR: %s\n',safeExceptionText(ME));
    rethrow(ME);
end

%% ---------------------- NORMAL SHUTDOWN ----------------------------
try
    if ~isempty(ioObj)
        io64(ioObj, LPT1address, 0);
    end
catch
end

try
    if exist('stim','var')
        safelyStopVNS(stim);
    end
catch
end

try
    if audio_active && ~isempty(pahandle), PsychPortAudio('Close', pahandle); end
catch
end

if eyelink_active
    try
        Eyelink('StopRecording');
        pause(1);
        Eyelink('CloseFile');
        Eyelink('ReceiveFile', edfFile, localEdfPath);
        Eyelink('Shutdown');
        disp('EyeLink session ended and EDF retrieved.');
    catch
        warning('EyeLink shutdown/retrieval problem.');
    end
end

%% ========================== FUNCTIONS ==============================



function site = stimConditionToSite(condition)
    switch upper(condition)
        case 'ACTIVE'
            site = 'LEFT_CYMBA_CONCHAE';
        case 'SHAM'
            site = 'LEFT_EARLOBE';
        otherwise
            site = 'NONE';
    end
end

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


function [stimPlan,matchedGSIdx,qc] = makeSessionVNSPlan(testIsRare,nGDstim,nGSstim,sessionCondition)
% Create the 180-trial test-phase stimulation plan for one visit.
% NEW HARD BUFFER RULE:
%   Every stimulation-paired trial (GD or matched GS) is immediately followed
%   by an UNSTIMULATED GS. That N+1 GS is flagged as a post-stimulation buffer.
%   Therefore neither another stimulation-paired trial nor a pure-tone GD can
%   immediately follow stimulation. Stimulation waveform/timing is unchanged.

    testIsRare=logical(testIsRare(:)');
    N=numel(testIsRare); rareIdx=find(testIsRare); standardIdx=find(~testIsRare);
    if N~=180 || numel(rareIdx)~=36 || nGDstim~=18 || nGSstim~=18
        error('Expected 180 test trials, 36 GD, 18 stimulated GD and 18 stimulated GS.');
    end
    sessionCondition=upper(sessionCondition);
    if ~ismember(sessionCondition,{'ACTIVE','SHAM'}), error('Session condition must be ACTIVE or SHAM.'); end

    % A trial can be stimulated only when its N+1 trial exists and is GS.
    % IMPORTANT: filter indices < N BEFORE indexing N+1. MATLAB evaluates
    % both sides of element-wise '&', so expressions such as
    % rareIdx<N & ~testIsRare(rareIdx+1) can still index N+1 when trial N
    % belongs to rareIdx.
    gdCandidates=rareIdx(rareIdx<N);
    gdCandidates=gdCandidates(~testIsRare(gdCandidates+1));
    if numel(gdCandidates)<nGDstim
        error('Not enough GD candidates followed by GS to enforce the post-stimulation buffer rule.');
    end

    % Search for a well-spread 18/18 GD split and 18 matched GS controls.
    bestScore=Inf; bestGD=[]; bestGS=[]; nAttempts=20000;
    for a=1:nAttempts
        gd=sort(gdCandidates(randperm(numel(gdCandidates),nGDstim)));
        % Stimulation events cannot be adjacent; their N+1 positions are buffers.
        if any(diff(gd)==1), continue; end
        gdBuffers=gd+1;

        % Stimulated GS must itself be followed by GS, cannot be a GD buffer,
        % and cannot be adjacent to any stimulated GD or to another stimulated GS.
        gsCand=standardIdx(standardIdx<N);
        gsCand=gsCand(~testIsRare(gsCand+1));
        gsCand=setdiff(gsCand,gdBuffers,'stable');
        gsCand=gsCand(~ismember(gsCand-1,gd) & ~ismember(gsCand+1,gd));
        if numel(gsCand)<nGSstim, continue; end
        ord=gsCand(randperm(numel(gsCand))); sel=[];
        for k=1:numel(ord)
            ix=ord(k);
            if isempty(sel) || all(abs(sel-ix)>1)
                % Its buffer cannot itself be a stimulation event or GD buffer.
                if ~ismember(ix+1,gd) && ~ismember(ix+1,gdBuffers) && ~ismember(ix+1,sel)
                    sel(end+1)=ix; %#ok<AGROW>
                    if numel(sel)==nGSstim, break; end
                end
            end
        end
        if numel(sel)~=nGSstim, continue; end
        gs=sort(sel); stim=sort([gd gs]); buffers=stim+1;
        if any(diff(stim)==1) || any(testIsRare(buffers)) || any(ismember(buffers,stim)), continue; end

        % Favor similar serial distributions for stimulated GD and GS.
        score=sum(abs(sort(gs)-sort(gd)));
        midpoint=(N+1)/2;
        score=score+2*abs(sum(gs<=midpoint)-sum(gd<=midpoint));
        if score<bestScore, bestScore=score; bestGD=gd; bestGS=gs; end
    end
    if isempty(bestGD)
        error('Could not construct the stimulation plan with mandatory post-stimulation GS buffers.');
    end

    stimulatedGDIdx=bestGD; matchedGSIdx=bestGS;
    pureToneGDIdx=setdiff(rareIdx,stimulatedGDIdx,'stable');
    stimMask=false(1,N); stimMask(stimulatedGDIdx)=true; stimMask(matchedGSIdx)=true;
    postStimBufferMask=false(1,N); postStimBufferMask(find(stimMask)+1)=true;

    % Hard QC: every stimulated trial must have exactly one clean GS N+1 buffer.
    stimIdx=find(stimMask); bufferIdx=find(postStimBufferMask);
    if any(stimIdx>=N) || any(testIsRare(stimIdx+1)) || any(stimMask(stimIdx+1))
        error('Internal VNS-plan error: mandatory unstimulated-GS N+1 buffer rule failed.');
    end
    if any(ismember(pureToneGDIdx,stimIdx+1))
        error('Internal VNS-plan error: a pure-tone GD immediately follows stimulation.');
    end
    if numel(bufferIdx)~=numel(stimIdx)
        error('Internal VNS-plan error: expected one unique buffer per stimulated trial.');
    end

    stimPlan=repmat({'NONE'},1,N); stimPlan(stimMask)={sessionCondition};
    qc=struct();
    qc.nGDTotal=numel(rareIdx); qc.nGDStimulated=numel(stimulatedGDIdx);
    qc.nGDPureTone=numel(pureToneGDIdx); qc.nGS=numel(matchedGSIdx);
    qc.nStimulated=sum(stimMask); qc.serialMatchScore=bestScore;
    qc.rareIdx=rareIdx; qc.stimulatedGDIdx=stimulatedGDIdx;
    qc.pureToneGDIdx=pureToneGDIdx; qc.matchedGSIdx=matchedGSIdx;
    qc.stimulatedIdx=stimIdx; qc.minTrialSeparation=min(diff(stimIdx));
    qc.noConsecutive=~any(diff(stimIdx)==1); qc.sessionCondition=sessionCondition;
    qc.postStimBufferIdx=bufferIdx; qc.postStimBufferMask=postStimBufferMask;
    qc.nPostStimBuffers=numel(bufferIdx);
    qc.allBuffersAreUnstimulatedGS=all(~testIsRare(bufferIdx) & ~stimMask(bufferIdx));
end

function labels = makeTimeBalancedLabels(indices, nActive, nSham)
% Randomly search balanced labels and minimize ACTIVE-vs-SHAM serial-position
% differences. Counts remain exact.

    indices = indices(:)';
    n = numel(indices);
    if nActive+nSham ~= n
        error('Label counts do not match number of indices.');
    end

    bestLabels = [];
    bestScore = Inf;

    for a = 1:2000
        tmp = [repmat({'ACTIVE'},1,nActive), repmat({'SHAM'},1,nSham)];
        tmp = tmp(randperm(n));

        activePos = indices(strcmp(tmp,'ACTIVE'));
        shamPos   = indices(strcmp(tmp,'SHAM'));

        % Primary balance: mean serial position.
        score = abs(mean(activePos)-mean(shamPos));

        % Secondary balance: counts in each half of the block.
        split = (min(indices)+max(indices))/2;
        score = score + 2*abs(sum(activePos<=split)-sum(shamPos<=split));

        if score < bestScore
            bestScore = score;
            bestLabels = tmp;
            if score == 0
                break;
            end
        end
    end

    labels = bestLabels;
end

function tone = makeChord(freqs, fs, toneDur, amplitude, rampDur)
% Create one 50-ms complex tone from three sinusoids.
% Rise/fall duration is 7 ms. The source paper specifies duration, not the
% mathematical envelope shape; this script retains the original linear ramp.

    nSamples = round(toneDur * fs);
    t = (0:nSamples-1) / fs;

    tone = zeros(1, nSamples);
    for f = freqs
        tone = tone + sin(2*pi*f*t);
    end

    tone = amplitude * tone / numel(freqs);

    rampSamples = round(rampDur * fs);
    if 2*rampSamples > nSamples
        error('Ramp duration is too long for the tone duration.');
    end

    rise = linspace(0, 1, rampSamples);
    fall = fliplr(rise);

    tone(1:rampSamples) = tone(1:rampSamples) .* rise;
    tone(end-rampSamples+1:end) = ...
        tone(end-rampSamples+1:end) .* fall;
end

function waveform = makePatternWaveform(pattern, toneA, toneB, fs, soa, toneDur)
% Build the complete 650-ms five-tone pattern as a single waveform.

    nPatternSamples = round((4*soa + toneDur) * fs);
    waveform = zeros(1, nPatternSamples);

    for s = 1:5
        startIdx = round((s-1) * soa * fs) + 1;

        if pattern(s) == 'A'
            tone = toneA;
        elseif pattern(s) == 'B'
            tone = toneB;
        else
            error('Pattern must contain only A and B.');
        end

        stopIdx = startIdx + length(tone) - 1;
        waveform(startIdx:stopIdx) = waveform(startIdx:stopIdx) + tone;
    end
end

function [isRare, preRareGaps, trailingStandards] = ...
    makeConstrainedRareSequence(nFrequent, nRare, minGap, maxGap)
% Generate a test sequence with:
%   - exactly nRare rare patterns,
%   - exactly nFrequent frequent patterns,
%   - minGap..maxGap standards before every rare.
%
% Therefore:
%   * the first rare follows 2-5 test standards,
%   * every later rare follows 2-5 standards after the preceding rare.
%
% Remaining standards, if any, occur after the final rare. Here we choose
% the trailing run randomly from 0..5 whenever feasible, keeping the end of
% the block from becoming unnecessarily long.

    minNeeded = nRare * minGap;
    maxPossible = nRare * maxGap;

    if nFrequent < minNeeded
        error('Not enough frequent trials to satisfy the minimum rare spacing.');
    end

    % Prefer 0-5 trailing standards, provided the pre-rare gaps can absorb
    % the remaining frequent trials.
    feasibleTail = [];
    for tail = 0:5
        targetGapSum = nFrequent - tail;
        if targetGapSum >= minNeeded && targetGapSum <= maxPossible
            feasibleTail(end+1) = tail; %#ok<AGROW>
        end
    end

    if isempty(feasibleTail)
        % General fallback if parameters are changed in the future
        minTail = max(0, nFrequent - maxPossible);
        maxTail = nFrequent - minNeeded;
        feasibleTail = minTail:maxTail;
    end

    trailingStandards = feasibleTail(randi(numel(feasibleTail)));
    targetGapSum = nFrequent - trailingStandards;

    % Start every gap at the minimum
    preRareGaps = minGap * ones(1, nRare);
    extra = targetGapSum - sum(preRareGaps);

    % Randomly distribute the remaining standards, never exceeding maxGap
    while extra > 0
        candidates = find(preRareGaps < maxGap);
        if isempty(candidates)
            error('Could not construct constrained rare sequence.');
        end

        idx = candidates(randi(numel(candidates)));
        preRareGaps(idx) = preRareGaps(idx) + 1;
        extra = extra - 1;
    end

    % Shuffle gap order so the composition itself is not systematic
    preRareGaps = preRareGaps(randperm(nRare));

    % Build logical sequence: false = frequent, true = rare
    isRare = false(1, nFrequent + nRare);
    pos = 1;

    for r = 1:nRare
        pos = pos + preRareGaps(r);
        isRare(pos) = true;
        pos = pos + 1;
    end

    % Validate counts and spacing
    if sum(isRare) ~= nRare
        error('Rare sequence has an incorrect number of rare trials.');
    end

    if sum(~isRare) ~= nFrequent
        error('Rare sequence has an incorrect number of frequent trials.');
    end

    rareIdx = find(isRare);

    % Number of standards before first rare
    firstGap = rareIdx(1) - 1;
    if firstGap < minGap || firstGap > maxGap
        error('First rare does not satisfy spacing constraint.');
    end

    % Standards between successive rares
    if numel(rareIdx) > 1
        between = diff(rareIdx) - 1;
        if any(between < minGap) || any(between > maxGap)
            error('Rare-spacing constraint failed.');
        end
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
