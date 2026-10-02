ER-taVNS ANALYSIS | README + MARKER MAP
=======================================
Use this as the single analysis note for the current 36-stimulated-GD +
36-pure-tone-GD artifact-buffer pilot.

RUN ORDER FOR FIRST PILOT
-------------------------
1. P0_ERtaVNS_marker_QC.m on E1 and E2 BrainVision .vmrk files.
2. B2_EEG_preparation_LG.m.
3. B3_EEG_preprocessing_LG.m.
4. Convert EyeLink EDF -> ASC, preserving the acquisition stem.
5. C1_Pupil_preprocessing.m.
6. Inspect marker QC, excluded-buffer counts, ordinary rejection counts and
   visit-specific plots before proceeding to group inference.

RAW MARKERS
-----------
S1 FLS; S2 FLD; S3 RLS; S4 RLD; S5 T1; S6 T2; S7 T3; S8 T4;
S9 T5 PRIMARY ANCHOR; S10 paired stimulation; E2: S11 circle, S12 square,
S13 triangle, S14 diamond, S15 space response.

POST-STIM BUFFER EEG ENCODING
-----------------------------
No dedicated EEG buffer TTL is sent. The buffer is defined deterministically:
trial N carrying S10 -> immediately following trial N+1 must be an unstimulated GS.
label_LG_events.m verifies this rule and marks N+1 as lg_post_stim_buffer.
EyeLink retains an explicit POST_STIM_BUFFER message as independent QC metadata.

NEW PRIMARY EXCLUSION RULE
--------------------------
Every stimulation-paired trial is followed by an unstimulated GS buffer.
Those N+1 buffer epochs are excluded from BOTH primary EEG and pupil analyses.
Raw data are never deleted. The flag remains available as:
  EEG: event/epoch lg_post_stim_buffer, lg_exclude_primary
  Pupil table: is_post_stim_buffer

EXPECTED COUNTS PER E1 OR E2 / VISIT
------------------------------------
410 T5 anchors (2 x 205 patterns); 360 TEST + 50 HABITUATION.
FLS/FLD/RLS/RLD = 169/169/36/36.
72 paired stimulation trials = 36 matched GS + 36 stimulated GD.
36 pure-tone GD. 72 post-stimulation GS buffers. 0 habituation stimulation.
After deterministic buffer exclusion, before ordinary artifact rejection:
TEST composites are expected to contain approximately LS=144, LD=144,
GS=216, GD=72. Exact later counts depend on ordinary blink/EEG rejection.

P0 MARKER QC
------------
Checks 410 tone markers per position, condition counts, 72 S10 trials,
72 N+1-derived buffers, no consecutive stimulation, and that every
stimulated trial is followed immediately by an unstimulated-GS buffer.

EEG PIPELINE
------------
label_LG_events.m: S9 remains the exact T5 anchor. It decodes S1..S4 and optional
S10, then derives N+1 buffers from the stimulation sequence without shifting event latencies.
set_LG_epoch_labels.m: transfers the zero-latency anchor metadata to each epoch.
B3_EEG_preprocessing_LG.m: excludes post-stim buffer epochs before ordinary EEG
artifact rejection, then continues the existing filtering/ICA/interpolation/
reference workflow. The stimulation-generation code is acquisition-side and was
not changed by this analysis update.

PUPIL PIPELINE
--------------
parse_LG_eyelink.m expects 410 trials and reads POST_STIM_BUFFER messages.
preprocess_pupil_LG.m excludes buffer trials from the primary clean table before
ordinary blink/missing-data rejection. C1 saves the original Local/Global plots
plus *_ERtaVNS_visit_overlay.png, which overlays the current visit's paired
(ACTIVE or SHAM) response. Visit 1 contains only whichever stimulation condition
was assigned to Visit 1; ACTIVE-vs-SHAM inference requires both visits.

E2 VISUAL/RESPONSE EVENTS
-------------------------
Acquisition sends S11/S12/S13/S14 at each circle/square/triangle/diamond onset
and S15 at each space response. EyeLink mirrors every visual onset and response,
including correct/false response metadata and RT. Auditory S5..S10 markers are
sent independently and remain T5-locked for Local/Global analysis.

LATER GROUP SCRIPTS
-------------------
B4-B7, C2-C3 and D* include legacy code from the previous study. They are useful
as templates/exports but are not by themselves the final confirmatory two-visit
Stimulation x Task x Global crossover model. Validate the pilot first.

VISIT-SPECIFIC PUPIL FIGURE (2026-09-21)
----------------------------------------
C1_Pupil_preprocessing.m automatically checks Visit1 and Visit2. Each existing
visit is processed separately and session_condition determines ACTIVE/taVNS vs
SHAM. The added Local and Global panels each contain four traces:
  Pure Standard, Pure Deviant, visit-condition + Standard, visit-condition + Deviant.
Pure = retained unstimulated TEST trials. Post-stimulation buffer trials are removed
before plotting. No cross-visit pooling is performed in this version.
