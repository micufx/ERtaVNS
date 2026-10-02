clc; clearvars; close all;

%% VOLUME CALIBRATION
% Participant-specific auditory level calibration.
%
% Run simply:
%       adjust_volume_LG
%
% No WAV file is required.
%
% CALIBRATION REFERENCE
% ---------------------
% One fixed low-pitched Local-Global chord is repeated throughout:
%
%       A = [350 700 1400] Hz --> Taken from Local-Global auditory paradigm
%
% Waveform parameters match the current Local-Global experiment:
%   sampling rate = 44100 Hz
%   tone duration = 50 ms
%   rise/fall     = 7 ms
%   amplitude     = 0.5
%   tone segment  = 50 ms tone + 100 ms silence
%
% PROCEDURE
% ---------
% 1. Set Windows master volume to 100.
% 2. Use the same audio output/headphones as in the experiment.
% 3. Keep the MATLAB "Volume calibration controls" window selected.
% 4. Adjust the repeated tone until it is JUST BARELY AUDIBLE.
% 5. Three refinement passes use the department's ORIGINAL linear gain steps:
%
%       0.025  ->  0.0025  ->  0.00025
%
% 6. For the final experiment-level calculation, the threshold gain is not
%    allowed to be lower than 0.00025:
%
%       thresholdGainForCalculation = max(measuredThreshold, 0.00025)
%
%    This preserves the practical lower resolution of the historical
%    calibration procedure and prevents unusually low subjective thresholds
%    from driving the experimental presentation level lower still.
%
% 7. The final Windows setting targets 40 dB sensation level (dB SL), i.e.
%    40 dB above the threshold used for calculation.
%
% KEYBOARD CONTROLS
% -----------------
%   LEFT ARROW  = quieter
%   RIGHT ARROW = louder
%   SPACE       = accept current level
%   ESC         = abort
%
% NOTE ABOUT PSYCHTOOLBOX
% -----------------------
% This version DOES NOT use Psychtoolbox KbCheck/KbWait/Screen().
% Keyboard input is handled with a normal MATLAB figure callback because
% Screen.mexw64 is not required for this calibration. PsychPortAudio is
% still used for sound playback.
%
% HISTORICAL WINDOWS-VOLUME TRANSFER
% ----------------------------------
% The department's fitted Windows-volume relationship is:
%
%       relative level = 0.9236 + 15.01 * ln(Windows slider)
%
% This transfer function is valid only for the audio-output setup for which
% it was measured.

% Modified version (16.09.2026) from Sebastian Puschmann (27.11.2024)

%% ---------------- Auditory parameters ---------------------
fs          = 44100;
toneDur     = 0.050;       % 50 ms
soa         = 0.150;       % 150 ms onset-to-onset
rampDur     = 0.007;       % 7 ms rise/fall
amplitude   = 0.5;         % same as experiment

freqs_A = [350 700 1400];  % taken from Local-Global auditory paradigm

referenceName  = 'A (low-pitched)';
referenceFreqs = freqs_A;

%% ---------------- Historical threshold procedure ----------------------
TARGET_SL_DB = 40;

gain_start    = 0.5;

% ORIGINAL department gain increments:
gain_stepsize = [0.025 0.0025 0.00025];

% Explicit lower bound used ONLY for the final threshold calculation.
% A participant may still adjust lower numerically if the interface allows
% it, but the calculated experiment level will never use a value < 0.00025.
thresholdFloor = 0.00025;

repeatPeriod = 1.0;        % approximately one presentation per second

% Historical Windows-volume fit
fit_a = 0.9236;
fit_b = 15.01;

% [] = current/default Windows output device
audioDevice = [];

%% ---------------- Generate reference stimulus --------------------------
tone = makeLGChord(referenceFreqs, fs, toneDur, amplitude, rampDur);

% 50 ms tone + 100 ms silence
silenceSamples = round((soa - toneDur) * fs);
referenceStim = [tone zeros(1, silenceSamples)];

peakSample = max(abs(referenceStim));
if peakSample >= 1
    error('Generated calibration waveform clips (peak = %.3f).', peakSample);
end

%% ---------------- Initialize PsychPortAudio ----------------------------
if exist('InitializePsychSound','file') == 0 || exist('PsychPortAudio','file') == 0
    error(['Psychtoolbox/PsychPortAudio was not found on the MATLAB path. ' ...
           'The calibration cannot run without PsychPortAudio.']);
end

InitializePsychSound(1);

% Close stale handles from previous failed attempts.
try
    PsychPortAudio('Close');
catch
end

% Match the experiment's basic audio configuration:
% playback mode, latency class 1, 44.1 kHz, stereo.
pahandle = PsychPortAudio('Open', audioDevice, 1, 1, fs, 2);
audioCleanup = onCleanup(@() closeLGCalibrationAudio(pahandle)); %#ok<NASGU>
PsychPortAudio('Volume', pahandle, 1);

%% ---------------- MATLAB keyboard-control window -----------------------
ctrlFig = figure( ...
    'Name','Volume calibration controls', ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Color',[0.94 0.94 0.94], ...
    'Position',[300 300 620 260], ...
    'WindowKeyPressFcn',@calibrationKeyPress);

figCleanup = onCleanup(@() closeControlFigure(ctrlFig)); %#ok<NASGU>

setappdata(ctrlFig,'LastKey','');

uicontrol(ctrlFig, ...
    'Style','text', ...
    'Units','normalized', ...
    'Position',[0.05 0.55 0.90 0.38], ...
    'BackgroundColor',[0.94 0.94 0.94], ...
    'FontSize',14, ...
    'HorizontalAlignment','center', ...
    'String',{ ...
        'VOLUME CALIBRATION', ...
        'Keep this window selected during calibration.', ...
        'LEFT = quieter     RIGHT = louder', ...
        'SPACE = accept     ESC = abort'});

statusTextHandle = uicontrol(ctrlFig, ...
    'Style','text', ...
    'Units','normalized', ...
    'Position',[0.05 0.15 0.90 0.25], ...
    'BackgroundColor',[0.94 0.94 0.94], ...
    'FontSize',13, ...
    'HorizontalAlignment','center', ...
    'String','Press any key in this window to start.');

figure(ctrlFig);
drawnow;

%% ---------------- Instructions ----------------------------------------
clc;
fprintf('============================================================\n');
fprintf(' VOLUME CALIBRATION\n');
fprintf('============================================================\n\n');
fprintf('Fixed reference tone: %s = [%s] Hz\n\n', ...
    referenceName, num2str(referenceFreqs));
fprintf('The SAME low-pitched reference tone is repeated throughout.\n');
fprintf('No WAV file is required.\n\n');
fprintf('1) Set the WINDOWS VOLUME to 100.\n');
fprintf('2) Use the same audio output/headphones as in the experiment.\n');
fprintf('3) Click once on the MATLAB window "Volume calibration controls".\n');
fprintf('4) Adjust the tone until it is JUST BARELY AUDIBLE.\n\n');
fprintf('LEFT ARROW  = quieter\n');
fprintf('RIGHT ARROW = louder\n');
fprintf('SPACE        = accept current level\n');
fprintf('ESC          = abort\n\n');
fprintf('Press any key IN THE CONTROL WINDOW to start.\n');

waitForAnyMatlabKey(ctrlFig);

%% ---------------- Three-step threshold refinement ----------------------
thres = nan(1,3);

for runIdx = 1:3
    if ~ishandle(ctrlFig)
        error('The control window was closed.');
    end

    if runIdx == 1
        volume = gain_start;
    else
        volume = thres(runIdx-1);
    end

    set(statusTextHandle,'String',sprintf([ ...
        'Pass %d/3   Gain step %.5f\n' ...
        'Press any key to begin this pass.'], ...
        runIdx, gain_stepsize(runIdx)));

    setappdata(ctrlFig,'LastKey','');
    figure(ctrlFig);
    drawnow;

    clc;
    fprintf('============================================================\n');
    fprintf(' PASS %d OF 3\n',runIdx);
    fprintf('============================================================\n\n');
    fprintf('Reference tone: %s = [%s] Hz\n\n', ...
        referenceName,num2str(referenceFreqs));
    fprintf('The same tone is played repeatedly.\n');
    fprintf('Adjust it until it is JUST BARELY AUDIBLE.\n\n');
    fprintf('LEFT ARROW  = quieter\n');
    fprintf('RIGHT ARROW = louder\n');
    fprintf('SPACE        = accept level\n');
    fprintf('ESC          = abort\n');
    fprintf('Gain step = %.5f\n\n',gain_stepsize(runIdx));
    fprintf('Press any key IN THE CONTROL WINDOW.\n');

    waitForAnyMatlabKey(ctrlFig);

    pause(0.4);
    done = false;

    while ~done
        if ~ishandle(ctrlFig)
            error('The control window was closed.');
        end

        % Preserve historical behavior: gain cannot go below zero.
        volume = min(max(volume,0),1);

        stereoStim = volume .* [referenceStim; referenceStim];
        PsychPortAudio('FillBuffer',pahandle,stereoStim);

        % Play one copy of the same fixed reference stimulus.
        PsychPortAudio('Start',pahandle,1,0,1);

        if volume > 0
            atten_dB = 20*log10(volume);
            atten_txt = sprintf('%.1f dB relative to experiment',atten_dB);
        else
            atten_txt = '-Inf dB (silence)';
        end

        set(statusTextHandle,'String',sprintf([ ...
            'Pass %d/3\n' ...
            'Gain = %.7f  (%s)\n' ...
            'LEFT quieter | RIGHT louder | SPACE accept'], ...
            runIdx,volume,atten_txt));

        fprintf('\rReference tone A | Gain = %.7f        ',volume);

        setappdata(ctrlFig,'LastKey','');
        t0 = tic;

        while ~done && toc(t0) <= repeatPeriod
            drawnow;

            if ~ishandle(ctrlFig)
                error('The control window was closed.');
            end

            key = getappdata(ctrlFig,'LastKey');

            if ~isempty(key)
                setappdata(ctrlFig,'LastKey','');

                switch lower(key)
                    case 'escape'
                        fprintf('\n');
                        error('Calibration was aborted with ESC.');

                    case 'leftarrow'
                        volume = max(0, volume - gain_stepsize(runIdx));
                        break;

                    case 'rightarrow'
                        volume = min(1, volume + gain_stepsize(runIdx));
                        break;

                    case 'space'
                        if volume <= 0
                            set(statusTextHandle,'String', ...
                                'Gain 0 is silence. Increase the level first.');
                            fprintf('\nGain 0 is silence. Increase the level first.\n');
                        else
                            done = true;
                        end
                end
            end

            pause(0.005);
        end
    end

    thres(runIdx) = volume;

    fprintf('\nPass %d saved: Gain = %.7f\n',runIdx,volume);
    set(statusTextHandle,'String',sprintf( ...
        'Pass %d saved: Gain = %.7f',runIdx,volume));

    pause(0.8);
end

%% ---------------- Threshold used for calculation -----------------------
measuredThreshold = thres(3);

% Preserve the practical lower resolution of the historical procedure.
thresholdGainForCalculation = max(measuredThreshold, thresholdFloor);

floorApplied = measuredThreshold < thresholdFloor;

%% ---------------- Convert threshold to Windows slider ------------------
% The normal experiment waveform is gain = 1 relative to the calibration
% waveform (both use amplitude = 0.5). Therefore:
SL_at_Windows100_dB = -20*log10(thresholdGainForCalculation);

% Required Windows-output change relative to slider 100.
requiredWindowsChange_dB = TARGET_SL_DB - SL_at_Windows100_dB;

% Historical Windows-slider relationship:
%       L(x) = fit_a + fit_b*ln(x)
fitLevelAt100  = fit_a + fit_b*log(100);
targetFitLevel = fitLevelAt100 + requiredWindowsChange_dB;
windows_slider_exact = exp((targetFitLevel-fit_a)/fit_b);

if windows_slider_exact > 100
    windows_slider_recommended = 100;
    statusResult = 'TARGET ABOVE AVAILABLE OUTPUT';
elseif windows_slider_exact < 1
    windows_slider_recommended = 1;
    statusResult = 'TARGET BELOW CALIBRATED WINDOWS RANGE';
else
    windows_slider_recommended = round(windows_slider_exact);
    statusResult = 'OK';
end

windowsChangeRounded_dB = fit_b*log(windows_slider_recommended/100);
predictedSL_rounded_dB = SL_at_Windows100_dB + windowsChangeRounded_dB;

%% ---------------- Results ----------------------------------------------
clc;
fprintf('============================================================\n');
fprintf(' AUDITORY CALIBRATION COMPLETE\n');
fprintf('============================================================\n\n');

fprintf('Calibration reference:                   %s [%s] Hz\n', ...
    referenceName,num2str(referenceFreqs));
fprintf('Measured final threshold gain:           %.7f\n',measuredThreshold);
fprintf('Minimum threshold used for calculation:  %.7f\n',thresholdFloor);
fprintf('Threshold gain used for calculation:     %.7f\n', ...
    thresholdGainForCalculation);

if floorApplied
    fprintf('Threshold floor applied:                 YES\n');
else
    fprintf('Threshold floor applied:                 NO\n');
end

fprintf('Experiment at Windows 100 would be:      %.2f dB above threshold\n', ...
    SL_at_Windows100_dB);
fprintf('Target experimental level:               %.1f dB SL\n',TARGET_SL_DB);
fprintf('Calculated Windows setting (exact):      %.2f\n',windows_slider_exact);
fprintf('Recommended Windows setting:             %d\n',windows_slider_recommended);
fprintf('Expected level after rounding:           %.2f dB SL\n', ...
    predictedSL_rounded_dB);
fprintf('Status:                                   %s\n\n',statusResult);

if floorApplied
    fprintf(['NOTE: The measured threshold was below 0.00025.\n' ...
             'For the final volume calculation, 0.00025 was used instead\n' ...
             'to preserve the lower resolution of the historical procedure.\n\n']);
end

if strcmp(statusResult,'OK')
    fprintf('=> Set the Windows volume to %d for the experiment.\n\n', ...
        windows_slider_recommended);
elseif windows_slider_exact > 100
    fprintf(['WARNING: The target of %.1f dB SL cannot be reached with the\n' ...
             'current audio chain at Windows volume 100.\n\n'],TARGET_SL_DB);
else
    fprintf(['WARNING: The calculated value is below the range of the historical\n' ...
             'Windows slider calibration.\n\n']);
end

fprintf(['IMPORTANT: Keep amplitude = 0.5 unchanged in the experiment and\n' ...
         'use the same audio output / headphones.\n']);

if ishandle(ctrlFig)
    set(statusTextHandle,'String',sprintf([ ...
        'Calibration finished.\n' ...
        'Recommended Windows volume = %d'], ...
        windows_slider_recommended));
end

%% ---------------- Leave results in workspace ---------------------------
LG_volume_calibration = struct();
LG_volume_calibration.referenceChord = 'A';
LG_volume_calibration.referenceDescription = 'low-pitched physical chord';
LG_volume_calibration.referenceFreqsHz = referenceFreqs;
LG_volume_calibration.thresholdRuns = thres;
LG_volume_calibration.measuredThresholdGain = measuredThreshold;
LG_volume_calibration.thresholdFloor = thresholdFloor;
LG_volume_calibration.thresholdFloorApplied = floorApplied;
LG_volume_calibration.thresholdGainUsedForCalculation = thresholdGainForCalculation;
LG_volume_calibration.adjustmentMethod = 'historical linear gain steps';
LG_volume_calibration.adjustmentSteps = gain_stepsize;
LG_volume_calibration.targetSL_dB = TARGET_SL_DB;
LG_volume_calibration.SL_at_Windows100_dB = SL_at_Windows100_dB;
LG_volume_calibration.windowsSliderExact = windows_slider_exact;
LG_volume_calibration.windowsSliderRecommended = windows_slider_recommended;
LG_volume_calibration.predictedSLAfterRounding_dB = predictedSL_rounded_dB;
LG_volume_calibration.status = statusResult;
LG_volume_calibration.fs = fs;
LG_volume_calibration.toneDur = toneDur;
LG_volume_calibration.soa = soa;
LG_volume_calibration.rampDur = rampDur;
LG_volume_calibration.amplitude = amplitude;
LG_volume_calibration.timestamp = datetime('now');

%% ======================= LOCAL FUNCTIONS ===============================

function tone = makeLGChord(freqs,fs,toneDur,amplitude,rampDur)
% Generate a Local-Global complex tone exactly as in the experiment.

    nSamples = round(toneDur*fs);
    t = (0:nSamples-1)/fs;

    tone = zeros(1,nSamples);
    for f = freqs
        tone = tone + sin(2*pi*f*t);
    end

    tone = amplitude*tone/numel(freqs);

    rampSamples = round(rampDur*fs);
    if 2*rampSamples > nSamples
        error('Ramp duration is too long for the tone duration.');
    end

    rise = linspace(0,1,rampSamples);
    fall = fliplr(rise);

    tone(1:rampSamples) = tone(1:rampSamples).*rise;
    tone(end-rampSamples+1:end) = ...
        tone(end-rampSamples+1:end).*fall;
end

function calibrationKeyPress(src,event)
    if ishandle(src)
        setappdata(src,'LastKey',event.Key);
    end
end

function waitForAnyMatlabKey(fig)
    if ~ishandle(fig)
        error('The control window was closed.');
    end

    setappdata(fig,'LastKey','');
    figure(fig);
    drawnow;

    while ishandle(fig)
        drawnow;
        key = getappdata(fig,'LastKey');

        if ~isempty(key)
            if strcmpi(key,'escape')
                error('Calibration was aborted with ESC.');
            end
            setappdata(fig,'LastKey','');
            return;
        end

        pause(0.01);
    end

    error('The control window was closed.');
end

function closeLGCalibrationAudio(pahandle)
    try
        PsychPortAudio('Stop',pahandle);
    catch
    end

    try
        PsychPortAudio('Close',pahandle);
    catch
    end
end

function closeControlFigure(fig)
    try
        if ishandle(fig)
            close(fig);
        end
    catch
    end
end
