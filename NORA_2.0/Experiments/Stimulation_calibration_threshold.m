%% Stimulation_calibration_threshold
% Interactive participant-controlled ER-taVNS intensity calibration.
%
% PURPOSE
% -------
% The participant adjusts stimulation strength using LEFT/RIGHT arrows.
% Each adjustment delivers ONE short stimulation train at the new level.
% SPACE starts confirmation; a second consecutive SPACE accepts the
% selected level. LEFT/RIGHT remain available during confirmation and any
% adjustment restarts the two-SPACE confirmation sequence. ESC aborts.
%
% IMPORTANT HARDWARE ASSUMPTION
% -----------------------------
% This script is configured for the current NORA laboratory setup:
%   NI device:       Dev2
%   analog output:   AO0 -> DS5 INPUT
%   digital gate:    port0/line0 (used if available)
%   NI rate:         200 kHz
%   stimulation:     25 Hz, 250 us pulses, 1.5 s train
%   command range:   0 to 5 V
%
% DS5_MA_PER_V = 1.0 below corresponds to the DS5 setting 10 V : 10 mA,
% i.e. approximately 1 V command = 1 mA output. If the DS5 range is changed,
% DO NOT use the displayed mA values until DS5_MA_PER_V is updated.
%
% SAFETY / PROCEDURE
% ------------------
% - Use only under the approved study protocol and established lab setup.
% - Verify electrode placement/contact and DS5 settings before starting.
% - The participant should select a level that is clearly perceptible and
%   strong but comfortable, NEVER painful.
% - LEFT immediately selects a weaker level for the next train.
% - The script never exceeds MAX_COMMAND_V.
% - First SPACE asks for confirmation and emits no additional stimulation.
% - Second consecutive SPACE accepts the level.
% - LEFT/RIGHT during confirmation adjust the level, deliver one new train,
%   and require the two-SPACE confirmation sequence again.
%
% The final command amplitude is the value to enter in the NORA experiment's
% existing 0-5 stimulation-amplitude prompt.
%
% Generated for the NORA 2.0 ER-taVNS workflow, 21.09.2026.

clc;
clearvars;
close all;
P = ERtaVNS_parameters(); % central editable stimulation/timing configuration

%% ---------------- User-adjustable constants ----------------------------
DEVICE_ID       = P.deviceID;
AO_CHANNEL      = P.aoChannel;
GATE_LINE       = P.gateLine;

FS              = P.niFs;
STIM_FREQ_HZ    = P.stimFrequencyHz;
PULSE_WIDTH_SEC = P.pulseWidthSec;
TRAIN_SEC       = P.stimDurationSec;

AMP_STEP_V      = P.calibrationStep;      % with current 10V:10mA DS5 setting ~= 0.2 mA
MIN_COMMAND_V   = P.minCommand;
MAX_COMMAND_V   = P.maxCommand;
START_COMMAND_V = P.startCommand;

% Current laboratory DS5 setting discussed for NORA:
DS5_MA_PER_V    = P.ds5mAPerV;      % 10 V input range : 10 mA output range

% A short zero tail guarantees return to zero after every train.
POST_ZERO_SEC   = P.postZeroSec;

%% ---------------- Basic checks ----------------------------------------
if exist('daqreset','file') == 0
    error(['MATLAB Data Acquisition Toolbox / daqreset was not found. ' ...
           'Cannot run stimulation calibration.']);
end

if exist('rectpuls','file') == 0
    error(['rectpuls() was not found. Signal Processing Toolbox is required ' ...
           'for the NORA-derived pulse generation.']);
end

if exist('KbCheck','file') == 0 || exist('KbName','file') == 0
    error(['Psychtoolbox keyboard functions were not found. ' ...
           'Add Psychtoolbox to the MATLAB path before calibration.']);
end

fprintf('============================================================\n');
fprintf(' ER-taVNS STIMULATION CALIBRATION\n');
fprintf('============================================================\n');
fprintf('NI pathway:       %s/AO%d -> DS5\n',DEVICE_ID,AO_CHANNEL);
fprintf('Train:            %.3f s | %.1f Hz | %.1f us\n', ...
    TRAIN_SEC,STIM_FREQ_HZ,PULSE_WIDTH_SEC*1e6);
fprintf('Adjustment step:  %.1f command V (~%.1f mA at %.2f mA/V)\n', ...
    AMP_STEP_V,AMP_STEP_V*DS5_MA_PER_V,DS5_MA_PER_V);
fprintf('Allowed range:    %.1f to %.1f command V\n',MIN_COMMAND_V,MAX_COMMAND_V);
fprintf('\n');
fprintf('IMPORTANT: displayed mA assumes the DS5 is set to 10 V : 10 mA.\n');
fprintf('Verify the DS5 range and electrode/contact setup before continuing.\n\n');

answer = questdlg( ...
    sprintf(['Confirm before calibration:\n\n' ...
             'DS5 setting is 10 V : 10 mA\n' ...
             'Electrodes/contact have been checked\n' ...
             'Participant understands: strong but comfortable, never painful']), ...
    'NORA stimulation calibration - hardware check', ...
    'Confirmed','Cancel','Cancel');

if ~strcmp(answer,'Confirmed')
    fprintf('Calibration cancelled before NI output.\n');
    return;
end

%% ---------------- Probe NI pathway ------------------------------------
daqreset;
probe = [];
hasGate = false;

try
    probe = daq.createSession('ni');
    probe.addAnalogOutputChannel(DEVICE_ID,AO_CHANNEL,'Voltage');

    try
        probe.addDigitalChannel(DEVICE_ID,GATE_LINE,'OutputOnly');
        hasGate = true;
    catch MEgate
        warning('Digital gate line unavailable; calibration will use AO0 only. %s', ...
            oneLine(MEgate.message));
        hasGate = false;
    end

    set(probe,'Rate',FS);
    probe.IsContinuous = false;
    release(probe);
    probe = [];
catch ME
    try
        if ~isempty(probe), release(probe); end
    catch
    end
    error('Could not initialize %s/AO%d. %s',DEVICE_ID,AO_CHANNEL,ME.message);
end

fprintf('NI stimulation output available: %s/AO%d\n',DEVICE_ID,AO_CHANNEL);
if hasGate
    fprintf('Digital gate available: %s\n',GATE_LINE);
else
    fprintf('Digital gate unavailable: AO0-only mode\n');
end

%% ---------------- Participant control window ---------------------------
currentAmp = START_COMMAND_V;
lastDeliveredAmp = NaN;
confirmed = false;
confirmPending = false;
aborted = false;
busy = false;
command = '';

fig = figure( ...
    'Name','Stimulation calibration', ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Color',[0.97 0.97 0.97], ...
    'WindowStyle','normal', ...
    'WindowState','fullscreen');

cleanupFig = onCleanup(@()safeCloseFigure(fig)); %#ok<NASGU>

% Main instruction.
annotation(fig,'textbox',[0.08 0.84 0.84 0.11], ...
    'String',{'STIMULATION CALIBRATION', ...
              'Choose a stimulation that is clearly perceptible and strong, but comfortable — NEVER painful.'}, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'EdgeColor','none');

% Sensation scale.
ax = axes('Parent',fig,'Position',[0.10 0.48 0.80 0.25]);
hold(ax,'on');
plot(ax,[0 100],[0 0],'-','LineWidth',4);
for k = 0:10:100
    plot(ax,[k k],[-0.08 0.08],'-','LineWidth',1.5);
    text(ax,k,-0.18,sprintf('%d',k),'HorizontalAlignment','center','FontSize',11);
end
plot(ax,[80 80],[-0.13 0.13],'--','LineWidth',2);
text(ax,80,0.22,'~80: strong but comfortable','HorizontalAlignment','center', ...
    'FontSize',13,'FontWeight','bold');
text(ax,0,0.22,'No sensation','HorizontalAlignment','left','FontSize',12);
text(ax,100,0.22,'Painful / too strong','HorizontalAlignment','right','FontSize',12);
xlim(ax,[-3 103]);
ylim(ax,[-0.35 0.38]);
axis(ax,'off');

instructionHandle = uicontrol(fig, ...
    'Style','text', ...
    'Units','normalized', ...
    'Position',[0.08 0.22 0.84 0.18], ...
    'BackgroundColor',[0.97 0.97 0.97], ...
    'FontSize',17, ...
    'HorizontalAlignment','center', ...
    'String',{ ...
        'RIGHT ARROW  = stronger (+0.2 mA)', ...
        'LEFT ARROW   = weaker (-0.2 mA)', ...
        'SPACE = review/confirm this level     ESC = stop calibration'});

statusHandle = uicontrol(fig, ...
    'Style','text', ...
    'Units','normalized', ...
    'Position',[0.08 0.08 0.84 0.10], ...
    'BackgroundColor',[0.97 0.97 0.97], ...
    'FontSize',15, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','center', ...
    'String','READY: Press RIGHT ARROW to deliver the first 0.2 mA stimulation train.');

figure(fig);
drawnow;

fprintf('\nParticipant controls:\n');
fprintf('  RIGHT = +%.1f (~%.1f mA)\n',AMP_STEP_V,AMP_STEP_V*DS5_MA_PER_V);
fprintf('  LEFT  = -%.1f (~%.1f mA)\n',AMP_STEP_V,AMP_STEP_V*DS5_MA_PER_V);
fprintf('  SPACE = confirmation (press SPACE again to accept)\n');
fprintf('  ESC   = abort\n\n');

%% ---------------- Interactive loop ------------------------------------
% Use Psychtoolbox keyboard polling rather than MATLAB figure callbacks.
% This is more reliable for a fullscreen participant display and is already
% available in the NORA experiment environment.
KbName('UnifyKeyNames');
keyRight = KbName('RightArrow');
keyLeft  = KbName('LeftArrow');
keySpace = KbName('space');
keyEsc   = KbName('ESCAPE');

% Ignore the key that may have been used to start the script/dialog.
KbReleaseWait;

while ishandle(fig) && ~confirmed && ~aborted
    drawnow;

    if busy
        WaitSecs(0.01);
        continue;
    end

    [keyIsDown,~,keyCode] = KbCheck;
    if ~keyIsDown
        WaitSecs(0.01);
        continue;
    end

    if keyCode(keyRight)
        thisCommand = 'right';
    elseif keyCode(keyLeft)
        thisCommand = 'left';
    elseif keyCode(keySpace)
        thisCommand = 'space';
    elseif keyCode(keyEsc)
        thisCommand = 'escape';
    else
        KbReleaseWait;
        continue;
    end

    % One physical key press = one calibration action. This prevents a held
    % arrow from stepping repeatedly or queuing multiple stimulation trains.
    KbReleaseWait;

    switch thisCommand
        case 'right'
            % Any intensity change cancels a pending confirmation. The
            % participant experiences the newly selected level and must then
            % press SPACE twice again to accept it.
            confirmPending = false;
            newAmp = min(MAX_COMMAND_V,currentAmp + AMP_STEP_V);
            if abs(newAmp-currentAmp) < 1e-12
                fprintf('[CAL] Maximum %.1f reached; no stronger train delivered.\n',MAX_COMMAND_V);
                set(statusHandle,'String',sprintf( ...
                    'Maximum reached: %.1f mA. SPACE starts confirmation; LEFT reduces.', ...
                    currentAmp*DS5_MA_PER_V));
                continue;
            end
            currentAmp = newAmp;

        case 'left'
            confirmPending = false;
            newAmp = max(MIN_COMMAND_V,currentAmp - AMP_STEP_V);
            if abs(newAmp-currentAmp) < 1e-12
                fprintf('[CAL] Minimum %.1f reached; no stimulation delivered.\n',MIN_COMMAND_V);
                set(statusHandle,'String','Minimum reached: 0.0 mA. RIGHT increases.');
                continue;
            end
            currentAmp = newAmp;

        case 'space'
            if ~confirmPending
                % First SPACE: ask for explicit confirmation, but do NOT
                % stimulate again. LEFT/RIGHT remain available and each
                % adjustment restarts the two-SPACE confirmation sequence.
                confirmPending = true;
                set(statusHandle,'String',sprintf( ...
                    ['Confirm %.1f mA?  SPACE = YES, FINAL CONFIRMATION   |   ' ...
                     'LEFT = weaker (-%.1f)   |   RIGHT = stronger (+%.1f)'], ...
                    currentAmp*DS5_MA_PER_V, ...
                    AMP_STEP_V*DS5_MA_PER_V,AMP_STEP_V*DS5_MA_PER_V));
                fprintf('[CAL] Confirm %.1f V (~%.1f mA): SPACE again = accept; LEFT/RIGHT = adjust.\n', ...
                    currentAmp,currentAmp*DS5_MA_PER_V);
                continue;
            else
                % Second consecutive SPACE: final acceptance.
                confirmed = true;
                break;
            end

        case 'escape'
            aborted = true;
            break;

        otherwise
            continue;
    end

    if confirmed || aborted
        break;
    end

    % Only LEFT/RIGHT reach this point. Deliver exactly one
    % experimental-style 1.5-s train at the newly selected level.
    busy = true;
    set(statusHandle,'String',sprintf( ...
        'Delivering %.1f mA ... please judge this level.', ...
        currentAmp*DS5_MA_PER_V));
    drawnow;

    fprintf('[CAL] Delivering command %.1f V (~%.1f mA) ...\n', ...
        currentAmp,currentAmp*DS5_MA_PER_V);

    try
        deliverCalibrationTrain(currentAmp,DEVICE_ID,AO_CHANNEL,GATE_LINE, ...
            hasGate,FS,STIM_FREQ_HZ,PULSE_WIDTH_SEC,TRAIN_SEC,POST_ZERO_SEC);
        lastDeliveredAmp = currentAmp;
        fprintf('[CAL] Train finished: %.1f V (~%.1f mA)\n', ...
            currentAmp,currentAmp*DS5_MA_PER_V);
    catch ME
        busy = false;
        set(statusHandle,'String','STIMULATION ERROR - calibration stopped.');
        drawnow;
        rethrow(ME);
    end

    busy = false;
    set(statusHandle,'String',sprintf( ...
        ['Current setting: %.1f mA.  LEFT = weaker | RIGHT = stronger | ' ...
         'SPACE = start confirmation'],currentAmp*DS5_MA_PER_V));
    figure(fig);
    drawnow;
end

%% ---------------- Finish ----------------------------------------------
if aborted || ~ishandle(fig)
    fprintf('\n============================================================\n');
    fprintf(' CALIBRATION ABORTED\n');
    fprintf('============================================================\n');
    fprintf('No threshold was accepted.\n');
    return;
end

% Require that a non-zero setting has actually been experienced.
if currentAmp <= 0 || isnan(lastDeliveredAmp)
    warning('No non-zero stimulation level was delivered and accepted.');
end

fprintf('\n============================================================\n');
fprintf(' CALIBRATION COMPLETE\n');
fprintf('============================================================\n');
fprintf('Selected command amplitude: %.1f\n',currentAmp);
fprintf('Approximate current at 10 V : 10 mA DS5 setting: %.1f mA\n', ...
    currentAmp*DS5_MA_PER_V);
fprintf('\nENTER THIS VALUE IN THE NORA EXPERIMENT:\n');
fprintf('    %.1f\n',currentAmp);
fprintf('============================================================\n');

set(statusHandle,'String',sprintf( ...
    'CONFIRMED: %.1f mA  |  NORA experiment amplitude = %.1f', ...
    currentAmp*DS5_MA_PER_V,currentAmp));
drawnow;

% Also leave convenient variables in the base workspace when run as script.
StimulationCalibration.commandAmplitude = currentAmp;
StimulationCalibration.approxCurrent_mA = currentAmp*DS5_MA_PER_V;
StimulationCalibration.ds5_mA_per_V = DS5_MA_PER_V;
StimulationCalibration.step = AMP_STEP_V;
StimulationCalibration.frequency_Hz = STIM_FREQ_HZ;
StimulationCalibration.pulseWidth_us = PULSE_WIDTH_SEC*1e6;
StimulationCalibration.train_s = TRAIN_SEC;
StimulationCalibration.timestamp = datetime('now');

assignin('base','StimulationCalibration',StimulationCalibration);

%% ======================== LOCAL FUNCTIONS ==============================
function deliverCalibrationTrain(amp,deviceID,aoChannel,gateLine,hasGate, ...
    fs,freqHz,pulseWidthSec,trainSec,postZeroSec)
% One finite NORA-derived train. The complete waveform is queued BEFORE
% output starts; no DataRequired callback is used.

    cycleTime = 1/freqHz;
    tCycle = -cycleTime/2:1/fs:cycleTime/2-(1/fs);
    cycleWave = amp * rectpuls(tCycle,pulseWidthSec);

    if amp ~= 0
        firstNZ = find(cycleWave ~= 0,1,'first');
        if isempty(firstNZ)
            error('Could not construct stimulation pulse.');
        end
        cycleWave = [cycleWave(firstNZ:end),cycleWave(1:firstNZ-1)];
    end

    nStimSamples = round(trainSec*fs);
    nCycles = ceil(nStimSamples/numel(cycleWave));
    train = repmat(cycleWave,1,nCycles);
    train = train(1:nStimSamples).';
    train(end) = 0;

    nPost = round(postZeroSec*fs);
    analogPacket = [train; zeros(nPost,1)];
    gatePacket = [ones(nStimSamples,1); zeros(nPost,1)];
    gatePacket(end) = 0;

    s = [];
    try
        s = daq.createSession('ni');
        s.addAnalogOutputChannel(deviceID,aoChannel,'Voltage');

        if hasGate
            s.addDigitalChannel(deviceID,gateLine,'OutputOnly');
            outputData = [analogPacket,gatePacket];
        else
            outputData = analogPacket;
        end

        set(s,'Rate',fs);
        s.IsContinuous = false;
        s.queueOutputData(outputData);

        % Foreground execution is intentional for calibration: one key press
        % corresponds to one complete train before another adjustment is read.
        s.startForeground();

        try
            release(s);
        catch
        end
    catch ME
        try
            if ~isempty(s)
                try, s.stop(); catch, end
                try, release(s); catch, end
            end
        catch
        end
        rethrow(ME);
    end
end

function safeCloseFigure(fig)
    try
        if ishandle(fig)
            close(fig);
        end
    catch
    end
end

function txt = oneLine(txt)
    txt = regexprep(char(txt),'[\r\n]+',' ');
end
