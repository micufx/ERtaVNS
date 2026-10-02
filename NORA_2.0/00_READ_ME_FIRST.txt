NORA 2.0 - CENTRAL PARAMETER CONFIGURATION
==========================================

The acquisition package now uses ONE central editable file:
    Experiments/ERtaVNS_parameters.m

For ordinary stimulation/timing changes, edit that file instead of A0/A1/B1.
A0 Control, A1 E1, B1 E2, and Stimulation_calibration_threshold all read it.

Example - test 1.0-s instead of 1.5-s stimulation:
    P.stimDurationSec = 1.000;

The calibration and all three acquisition scripts then use 1.0 s automatically.
The default supplied in this package remains the current 1.500-s protocol.

Minimum stimulation onset spacing is controlled independently by:
    P.minimumStimOnsetSec = 5.000;

Because every paired stimulation trial in E1/E2 has a mandatory N+1 GS buffer,
the earliest next stimulation is N+2. The script therefore derives the LG
silent-gap range automatically. With the present 650-ms pattern and a 5.0-s
minimum onset spacing, the derived gap remains 1.850-2.150 s (300-ms jitter),
regardless of whether the train itself is 1.0 or 1.5 s. This is mathematically
correct for an ONSET-to-ONSET requirement.

For Control, where there is no intervening LG buffer trial, the post-train OFF
time is derived from stimulation duration so consecutive stimulation onsets still
obey the same 5.0-s minimum plus 300-ms jitter.

The parameter function validates incompatible values and stops with an explicit
error instead of silently running an impossible timing combination.

IMPORTANT: Changing stimulation parameters changes the study protocol. Keep the
values consistent with the approved protocol and bench-check DS5 timing/output
before participant acquisition.

---------------------------------------------------------------------

NORA 2.0 — ULTIMATE READY (2026-09-21)
========================================
This package was rebuilt from NORA 2.0_PILOT_READY and contains both complete
Experiments and Analysis folders, including the original hardware/support files.

Main protocol changes retained:
- 25 habituation + 180 test trials per block; 2 blocks/task/visit.
- 36 GD per block: 18 stimulation-paired + 18 pure-tone GD.
- 18 matched GS stimulation-paired per block.
- Every stimulation-paired trial is immediately followed by an unstimulated GS
  buffer trial. The whole N+1 five-tone trial is excluded from primary EEG/pupil.
- Physical stimulation generation/delivery code is otherwise preserved.
- E2 markers: tones S5-S9; condition S1-S4; paired S10; figures S11-S14;
  response S15. Buffer is a repeated same S1-S4 condition code + EyeLink message.
- Participant IDs no longer abort simply because data already exist; operator can
  keep the ID or accept the suggested next unused participant.
- Unexpected MATLAB errors now display/rethrow the ORIGINAL error after emergency
  saving, rather than hiding it behind a generic message.
- Pupil analysis keeps T5 as the true epoch anchor, excludes buffer epochs, and
  automatically processes Visit1/Visit2 separately. Each visit gets a four-trace
  Local/Global figure: Pure Standard, Pure Deviant, and the visit's ACTIVE/taVNS
  Standard+Deviant or SHAM Standard+Deviant curves.

Important validation boundary:
The code and existing recorded marker files were inspected statically. MATLAB,
BrainVision, EyeLink, NI/DS5 and the lab parallel-port hardware are not available
in this environment, so a real hardware acquisition cannot be executed here.
Run a short lab test before participant acquisition and inspect the resulting .vmrk.
