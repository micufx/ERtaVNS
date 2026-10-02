function P = ERtaVNS_parameters()
% ERtaVNS_PARAMETERS - single editable parameter file for NORA 2.0.
%
% EDIT ONLY THE VALUES IN THE USER-EDITABLE SECTION unless you intentionally
% redesign the paradigm. A0/A1/B1 and Stimulation_calibration_threshold all
% read this file, so stimulation duration/frequency/pulse width stay matched.
%
% Timing rule used by A1/B1:
%   stimulation can occur no earlier than trial N+2 because N+1 is the
%   mandatory unstimulated GS buffer. Therefore, to guarantee a minimum
%   stimulation-ONSET-to-stimulation-ONSET interval M, the minimum pattern
%   onset interval is M/2. The silent gap is derived automatically.

%% ================= USER-EDITABLE PARAMETERS ======================
P.stimDurationSec       = 1.500;   % CURRENT protocol; set 1.000 to test a 1-s train
P.stimFrequencyHz       = 25;
P.pulseWidthSec         = 250e-6;
P.minimumStimOnsetSec   = 5.000;   % minimum onset-to-onset interval
P.ITIJitterSec          = 0.300;   % total silent-gap jitter range
P.ITIStepSec            = 0.050;   % discrete jitter step in E1/E2
P.calibrationStep       = 0.2;     % NORA 0-5 command convention
P.minCommand            = 0.0;
P.maxCommand            = 5.0;
P.startCommand          = 0.0;

% NI / DS5 setup (change only if hardware setup changes)
P.deviceID              = 'Dev2';
P.aoChannel             = 0;
P.gateLine              = 'port0/line0';
P.niFs                  = 200000;
P.ds5mAPerV             = 1.0;     % valid for current 10 V : 10 mA DS5 setting
P.postZeroSec           = 1.500;   % packet zero tail; does NOT change stim duration

% Local-global timing. Normally leave these unchanged.
P.audioFs               = 44100;
P.toneDurSec            = 0.050;
P.soaSec                = 0.150;
P.taskStartLeadSec      = 0.400;
%% ================================================================

% ---------------- Derived values: DO NOT hand-edit ----------------
P.patternDurSec = 4*P.soaSec + P.toneDurSec;
P.t5FromT1Sec = 4*P.soaSec;
% Packet must start taskStartLead before T1 and become non-zero at T5.
P.preZeroSec = P.taskStartLeadSec + P.t5FromT1Sec;

% Mandatory N+1 buffer means earliest next stimulation is N+2.
P.minSeriesOnsetSec = P.minimumStimOnsetSec/2;
P.silentGapMinSec = P.minSeriesOnsetSec - P.patternDurSec;
P.silentGapMaxSec = P.silentGapMinSec + P.ITIJitterSec;
P.silentGapValuesSec = P.silentGapMinSec:P.ITIStepSec:P.silentGapMaxSec;
P.seriesOnsetIntervalValuesSec = P.patternDurSec + P.silentGapValuesSec;

% Control has no intervening LG buffer trial, so derive its OFF time so
% consecutive control stimulation ONSETS obey the same minimum + jitter.
P.controlOffMinSec = P.minimumStimOnsetSec - P.stimDurationSec;
P.controlOffMaxSec = P.controlOffMinSec + P.ITIJitterSec;

% ---------------- Safety/consistency validation --------------------
assert(P.stimDurationSec > 0,'stimDurationSec must be > 0.');
assert(P.stimFrequencyHz > 0,'stimFrequencyHz must be > 0.');
assert(P.pulseWidthSec > 0 && P.pulseWidthSec < 1/P.stimFrequencyHz, ...
    'pulseWidthSec must be >0 and shorter than one stimulation cycle.');
assert(P.minimumStimOnsetSec > 0,'minimumStimOnsetSec must be > 0.');
assert(P.stimDurationSec < P.minimumStimOnsetSec, ...
    'Stimulation duration must be shorter than minimum stimulation onset interval.');
assert(P.silentGapMinSec >= 0, ...
    'Requested minimum stimulation spacing is incompatible with the LG pattern duration.');
assert(P.ITIJitterSec >= 0 && P.ITIStepSec > 0,'ITI jitter/step values are invalid.');
assert(P.maxCommand >= P.minCommand && P.calibrationStep > 0,'Calibration limits/step invalid.');

% Numerical check of the actual discrete timing grid.
P.actualMinimumStimOnsetSec = 2*min(P.seriesOnsetIntervalValuesSec);
if P.actualMinimumStimOnsetSec + 1e-9 < P.minimumStimOnsetSec
    error('Derived LG timing violates minimum stimulation onset interval.');
end
end
