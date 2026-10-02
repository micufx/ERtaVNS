clc; clear; close all;
% Updated C1: mean +/- SEM, buffer checks, and visit-matched control overlay.

%%
% Keep this corrected workflow folder first on the MATLAB path.
LG_THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
if ~isempty(LG_THIS_SCRIPT_DIR)
    addpath(LG_THIS_SCRIPT_DIR, '-begin');
end

%% Pupil preprocessing (batch) - Local/Global (OUTPUT epoch locked to TONE5)
% Corrected trigger workflow:
%   - EyeLink LG_TONE / LG_T5 messages are written at tone onset.
%   - No 50-ms subtraction is performed.
%   - T5 messages already contain the full Local-Global condition.
%
% OUTPUT epoch (saved + plotted):
%   - ALWAYS relative to TONE5 onset (time 0 = tone5 onset)
%   - ALWAYS fixed to epochWin_out_ms = [-500 2000]
%
% Baseline (ONLY thing selectable):
%   baselineAnchor = 'TONE5' or 'TONE1'
%   baselineWin_ms = e.g. [-500 0] relative to selected tone
%
% If baselineAnchor='TONE1', baseline may occur before -500 ms in TONE5 coords.
% We handle this INTERNALLY by extracting a bigger window only for baseline,
% then we crop back to [-500 2000] for output (so plots/grand-average stay consistent).
%
% Timing QC:
%   - SOA = diffs within trial (T1..T5)
%   - ITI = TONE1(next pattern) - TONE1(current pattern): 2500-2800 ms jitter
% Stores into T_LG:
%   SOA_E1_Pupil, SOA_E2_Pupil, ITI_E1_Pupil, ITI_E2_Pupil
%   BaselineMapped_E1_Pupil, BaselineMapped_E2_Pupil (string like "[-2160 -1660]")
%   EpochSuggest_E1_Pupil, EpochSuggest_E2_Pupil     (string like "[-2200 2000]")
%
% Miguel Contreras-Altamirano, 2026

%% Paths
PROJECT_ROOT = fileparts(LG_THIS_SCRIPT_DIR);
path = fullfile(PROJECT_ROOT,'Experiments','Data');

% AUTO VISIT MODE: process every Visit1 / Visit2 folder that exists.
% Each visit is kept separate. The visit condition (ACTIVE or SHAM) is read
% from the trial metadata parsed from the EyeLink file; it is NOT inferred
% from visit number.
visitsToAnalyze = [1 2];

%% ================= USER SETTINGS =================
startSubIdx = 1;           % set e.g. 10 to resume
endSubIdx   = inf;

fsPupil   = 500;           % Hz
dt        = 1000/fsPupil;  % ms (2 ms)
pad_ms    = 50;            % blink padding

% The corrected experiment writes EyeLink messages at tone onset, so no
% post-hoc 50-ms timestamp correction is needed.

% OUTPUT epoch fixed (for saving + plotting + later grand avg)
epochWin_out_ms = [-500 2000];  % form -500 in papers
relTimes_out    = epochWin_out_ms(1):dt:epochWin_out_ms(2);

% Filter
lpCut = 5;
[b_lp, a_lp] = butter(2, lpCut/(fsPupil/2), 'low');

blinkLimit = 0.20;

% Baseline selection (ONLY baseline changes)
baselineAnchor = 'TONE5';    % 'TONE5' or 'TONE1'
baselineWin_ms = [-500 0];   % relative to selected tone

% ================= PLOT / CONTROL SETTINGS =================
pupilUpperLimit = 1.2;      % Exact upper edge for EVERY pupil panel.
controlBaselineWin_ms = baselineWin_ms; % Same window, relative to stimulation onset.
% Control is preprocessed once per visit with the same blink/filter rules.
% For each overlay, scale control with THAT TASK's mean/SD used by PUPIL_Z;
% do not independently rescale the fifth curve to a different z-score unit.
assert(isscalar(pupilUpperLimit) && isfinite(pupilUpperLimit),'Invalid pupil upper limit.');

%% --- Find participants ---
dSubs   = dir(path);
isSub   = [dSubs.isdir] & ~ismember({dSubs.name},{'.','..'});
subNames = {dSubs(isSub).name};
subNames = subNames(~cellfun('isempty', regexp(subNames,'^S\d+','once')));

if isempty(subNames), error('No participant folders found under: %s', path); end
if isinf(endSubIdx), endSubIdx = numel(subNames); end
endSubIdx   = min(endSubIdx, numel(subNames));
startSubIdx = max(1, startSubIdx);

%% ================= MAIN LOOP =================
for analysisVisit = visitsToAnalyze

    outpath = fullfile(PROJECT_ROOT,'Analysis_Output',sprintf('Visit%d',analysisVisit));
    infoFile = fullfile(outpath,'Info_EEG_LG.mat');
    if ~exist(outpath,'dir'), mkdir(outpath); end
    haveTimingInfo = exist(infoFile,'file')==2;
    if haveTimingInfo
        load(infoFile,'T_LG');
        T_LG = ensure_pupil_timing_columns(T_LG);
    else
        fprintf('\n[Visit %d] No EEG info table; pupil/control preprocessing will still run.\n',analysisVisit);
    end

    fprintf('\n#################### VISIT %d ####################\n',analysisVisit);

    for sub = startSubIdx:endSubIdx

        participant = subNames{sub};
        subFolder   = fullfile(path,participant,sprintf('Visit%d',analysisVisit));
        if ~exist(subFolder,'dir')
            fprintf('   [%s] Visit%d folder absent. Skipping.\n',participant,analysisVisit);
            continue;
        end

        fprintf('\n====================================================\n');
        fprintf('Participant: %s (%d/%d) | Visit %d\n', participant, sub, numel(subNames), analysisVisit);
        fprintf('Folder:      %s\n', subFolder);
        fprintf('====================================================\n');

        outSubFolder = fullfile(outpath, participant);
        if ~exist(outSubFolder,'dir'), mkdir(outSubFolder); end

        % Same participant and visit only; never borrow another visit's control.
        control = struct('available',false,'condition',"",'dataTableClean',table, ...
            'source',"",'reason',"Control ASC not found");
        try
            dControl=dir(fullfile(subFolder,sprintf('%sV%dC*.asc',participant,analysisVisit)));
            exact=sprintf('%sV%dC.asc',participant,analysisVisit);
            match=find(strcmpi({dControl.name},exact));
            if ~isempty(match), dControl=dControl(match); end
            assert(numel(dControl)<=1,'Several control ASC files match; keep one exact normal stem.');
            if ~isempty(dControl)
                control.source=string(fullfile(subFolder,dControl(1).name));
                fprintf('   [%s V%d CONTROL] Using ASC: %s\n',participant,analysisVisit,control.source);
                controlLines=regexp(fileread(control.source),'\r\n|\n|\r','split');
                control=preprocess_control_pupil(controlLines,analysisVisit,control.source, ...
                    fsPupil,dt,pad_ms,relTimes_out,b_lp,a_lp,blinkLimit,controlBaselineWin_ms);
                controlDataTableClean=control.dataTableClean;
                controlRemovedTrials=control.removedTrials;
                controlQC=control.qc; relTimes=relTimes_out;
                save(fullfile(outSubFolder,sprintf('%s_V%d_C_pupil_preproc.mat',participant,analysisVisit)), ...
                    'controlDataTableClean','controlRemovedTrials','controlQC','relTimes','analysisVisit','-v7.3');
                fprintf('   [%s V%d CONTROL] %d/%d bursts retained (%s).\n',participant,analysisVisit, ...
                    height(controlDataTableClean),controlQC.n_bursts,control.condition);
            end
        catch ME
            control.available=false; control.reason=string(ME.message);
            warning('C1:ControlSkipped','%s V%d control: %s',participant,analysisVisit,ME.message);
        end
        if ~control.available
            fprintf('   [%s V%d] Control overlay unavailable: %s\n',participant,analysisVisit,control.reason);
        end

        for eidx = 1:2
            thisExp = sprintf('E%d', eidx);

            dAsc=dir(fullfile(subFolder,sprintf('%sV%d%s*.asc',participant,analysisVisit,thisExp)));
            dAsc=prefer_exact_stem(dAsc,sprintf('%sV%d%s.asc',participant,analysisVisit,thisExp));
            if isempty(dAsc)
                fprintf('   [%s V%d] No %sV%d%s*.asc found. Skipping.\n',participant,analysisVisit,participant,analysisVisit,thisExp);
                continue;
            end

            inputFile = fullfile(subFolder, dAsc(1).name);
            fprintf('   [%s | V%d | %s] Using ASC: %s\n', participant,analysisVisit,thisExp,dAsc(1).name);

            try
            raw   = fileread(inputFile);
            lines = regexp(raw, '\r\n|\n|\r', 'split');
            lg = parse_LG_eyelink(lines, 'expectedTrials', 410, 'verbose', true);
            tone = lg.tone;
            % Revalidate the COMPLETE raw trial sequence, before any rejection.
            lg=verify_buffer_sequence(lg);
            inputCondition=get_visit_condition(lg.trialInfo);
            plotControl=control;
            if plotControl.available && plotControl.condition~=inputCondition
                plotControl.available=false;
                plotControl.reason="Control ACTIVE/SHAM metadata do not match this task";
            end
            v=lg.trialInfo.visit(isfinite(lg.trialInfo.visit));
            assert(all(v==analysisVisit),'EyeLink visit metadata do not match the selected visit folder.');

            if haveTimingInfo

            [timingPupil, T_LG] = check_effective_timing_pupil( ...
                participant, T_LG, thisExp, tone, ...
                'baselineAnchor', baselineAnchor, ...
                'baselineWin_ms', baselineWin_ms, ...
                'epochEnd_ms', epochWin_out_ms(2), ...
                'store', 'median', ...
                'verbose', true);

            save(infoFile,'T_LG','-append');
            excelFile = fullfile(outpath,'Info_EEG_LG.xlsx');
            write_TLG_excel_exact(T_LG, excelFile);
            end
                results = preprocess_pupil_LG( ...
                    lines, lg, fsPupil, dt, pad_ms, relTimes_out, ...
                    b_lp, a_lp, blinkLimit, ...
                    'baselineAnchor', baselineAnchor, ...
                    'baselineWin_ms', baselineWin_ms);

                dataTableClean = results.dataTableClean;
                assert_buffer_exclusion(dataTableClean,lg.trialInfo);
                bufferQC=lg.report;
                bufferQC.n_retained=height(dataTableClean);
                bufferQC.n_retained_buffers=nnz(dataTableClean.is_post_stim_buffer);
                removedTrialsSummaryAll = results.removedTrialsSummaryAll;
                relTimes = relTimes_out;

                outMat = fullfile(outSubFolder,sprintf('%s_V%d_%s_pupil_preproc.mat',participant,analysisVisit,thisExp));
                save(outMat,'dataTableClean','removedTrialsSummaryAll','relTimes','analysisVisit','bufferQC','-v7.3');

                if isempty(dataTableClean) || height(dataTableClean)==0
                    warning('   [%s | V%d | %s] No retained trials. Skipping plots.',participant,analysisVisit,thisExp);
                    continue;
                end

                % Original descriptive figures are retained.
                fig1=plot_bekinschtein_pupil(dataTableClean,relTimes,'PUPIL_Z',participant,thisExp);
                apply_pupil_limits(fig1,pupilUpperLimit);
                saveas(fig1,fullfile(outSubFolder,sprintf('%s_V%d_%s_meanSEM_bySound.png',participant,analysisVisit,thisExp))); close(fig1);
                fig2=plot_mz4_pupil(dataTableClean,relTimes,'PUPIL_Z',participant,thisExp);
                apply_pupil_limits(fig2,pupilUpperLimit);
                saveas(fig2,fullfile(outSubFolder,sprintf('%s_V%d_%s_meanSEM_byType.png',participant,analysisVisit,thisExp))); close(fig2);
                fig3=plot_bekinschtein_dualpanel_pupil(dataTableClean,relTimes,'PUPIL_Z',participant,thisExp);
                apply_pupil_limits(fig3,pupilUpperLimit);
                saveas(fig3,fullfile(outSubFolder,sprintf('%s_V%d_%s_meanSEM_LocalGlobal.png',participant,analysisVisit,thisExp))); close(fig3);

                % Mean +/- SEM: four task curves plus the same-visit control.
                %   Pure Standard / Pure Deviant = unstimulated, non-buffer trials
                %   <ACTIVE|SHAM> + Standard / Deviant = stimulation-paired trials
                % The visit condition is read from session_condition in the data.
                fig4=plot_ERtaVNS_visit_fourtrace(dataTableClean,relTimes,'PUPIL_Z',participant,thisExp,analysisVisit,plotControl);
                controlOverlay=getappdata(fig4,'controlOverlay');
                save(outMat,'controlOverlay','-append');
                condLabel=get_visit_condition(dataTableClean);
                apply_pupil_limits(fig4,pupilUpperLimit);
                saveas(fig4,fullfile(outSubFolder,sprintf('%s_V%d_%s_Pure_plus_%s_LocalGlobal.png',participant,analysisVisit,thisExp,char(condLabel))));
                close(fig4);

                fprintf('   [%s | V%d | %s] Plots saved (%s).\n',participant,analysisVisit,thisExp,char(condLabel));

            catch ME
                fprintf(2,'   [%s | V%d | %s] ERROR: %s\n',participant,analysisVisit,thisExp,ME.message);
                disp(getReport(ME,'extended'));
                continue;
            end
        end
    end
end

fprintf('\nAll done.\n');
disp('Pupil/control processing finished; existing EEG timing tables updated when available.');


%% ==================== FUNCTIONS ========================

%% Figure Local-Global

function fig = plot_bekinschtein_dualpanel_pupil(T, relTimes, fieldName, participant, thisExp)
% plot_bekinschtein_dualpanel_pupil
% Creates a single figure with two panels:
%   TOP: Local Standard vs Local Deviant   (collapsed across global status)
%   BOT: Global Standard vs Global Deviant (collapsed across local status)
%
% Uses your existing labels in dataTableClean.type:
%   "Frequent Local Standard"
%   "Frequent Local Deviant"
%   "Rare Local Standard"
%   "Rare Local Deviant"
%
% Legend FIX:
%   - SEM patches are hidden from legend
%   - legend is created explicitly from line handles only
%
% Labels:
%   - No "(Rare/Frequent)" in legend text (your request)

typeVec = string(T.type);

% --- Define the collapsed factorial contrasts correctly ---
% Local Standard  = FLS + RLS
% Local Deviant   = FLD + RLD
% Global Standard = FLS + FLD
% Global Deviant  = RLS + RLD
idx_LS = (typeVec == "Frequent Local Standard") | ...
         (typeVec == "Rare Local Standard");
idx_LD = (typeVec == "Frequent Local Deviant") | ...
         (typeVec == "Rare Local Deviant");
idx_GS = (typeVec == "Frequent Local Standard") | ...
         (typeVec == "Frequent Local Deviant");
idx_GD = (typeVec == "Rare Local Standard") | ...
         (typeVec == "Rare Local Deviant");

% --- Helper to compute mean+SEM ---
    function [mu, sem, nTr] = mean_sem_from_idx(idx)
        nTr = sum(idx);
        if nTr == 0
            mu  = nan(1, numel(relTimes));
            sem = nan(1, numel(relTimes));
            return;
        end
        X = vertcat(T.(fieldName){idx}); % [nTrials x nTime]
        mu = mean(X, 1, 'omitnan');
        sd = std(X, 0, 1, 'omitnan');
        nT = sum(~isnan(X), 1);
        sem = sd ./ sqrt(max(nT, 1));
    end

[mu_LS, sem_LS, nLS] = mean_sem_from_idx(idx_LS);
[mu_LD, sem_LD, nLD] = mean_sem_from_idx(idx_LD);
[mu_GS, sem_GS, nGS] = mean_sem_from_idx(idx_GS);
[mu_GD, sem_GD, nGD] = mean_sem_from_idx(idx_GD);

% Match COMPARE figure aspect (wide + shorter height)
fig = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 0.65], ...
             'InvertHardcopy','off');

% ===================== TOP: LOCAL =====================
ax1 = subplot(1,2,1); hold(ax1,'on');
set(ax1,'Color','w','XColor','k','YColor','k');

title(ax1, sprintf('%s | %s | Local | Ntrials Std=%d Dev=%d', upper(string(participant)), upper(string(thisExp)), nLS, nLD), FontSize=12);
xlabel(ax1,'Time (ms)', FontSize=12);
ylabel(ax1,'Pupil size (z-score)', FontSize=12);

% SEM patches (hidden from legend)
p1 = fill(ax1, [relTimes fliplr(relTimes)], [(mu_LS+sem_LS) fliplr(mu_LS-sem_LS)], [0 0.4470 0.7410], ...
    'FaceAlpha',0.15,'EdgeColor','none');
set(p1,'HandleVisibility','off');

p2 = fill(ax1, [relTimes fliplr(relTimes)], [(mu_LD+sem_LD) fliplr(mu_LD-sem_LD)], [0.8500 0.3250 0.0980], ...
    'FaceAlpha',0.15,'EdgeColor','none');
set(p2,'HandleVisibility','off');

% Lines (these are the ONLY things used in legend)
hLS = plot(ax1, relTimes, mu_LS, 'LineWidth',2, 'Color',[0 0.45 0.85]);
hLD = plot(ax1, relTimes, mu_LD, 'LineWidth',2, 'Color',[0.90 0.45 0.10]);

xline(ax1, 0, '-', 'Color',[0.5 0.5 0.5]);
grid(ax1,'on');

leg1 = legend(ax1, [hLS hLD], {'Local Standard','Local Deviant'}, 'Location','southwest', FontSize=12);
set(leg1,'AutoUpdate','off');

hold(ax1,'off');

% ===================== BOTTOM: GLOBAL =====================
ax2 = subplot(1,2,2); hold(ax2,'on');
set(ax2,'Color','w','XColor','k','YColor','k');

title(ax2, sprintf('%s | %s | Global | Ntrials Std=%d Dev=%d', upper(string(participant)), upper(string(thisExp)), nGS, nGD), FontSize=12);
xlabel(ax2,'Time (ms)', FontSize=12);
ylabel(ax2,'Pupil size (z-score)', FontSize=12);

p3 = fill(ax2, [relTimes fliplr(relTimes)], [(mu_GS+sem_GS) fliplr(mu_GS-sem_GS)], [0 0.4470 0.7410], ...
    'FaceAlpha',0.15,'EdgeColor','none');
set(p3,'HandleVisibility','off');

p4 = fill(ax2, [relTimes fliplr(relTimes)], [(mu_GD+sem_GD) fliplr(mu_GD-sem_GD)], [0.8500 0.3250 0.0980], ...
    'FaceAlpha',0.15,'EdgeColor','none');
set(p4,'HandleVisibility','off');

hGS = plot(ax2, relTimes, mu_GS, 'LineWidth',2, 'Color',[0 0.45 0.85]);
hGD = plot(ax2, relTimes, mu_GD, 'LineWidth',2, 'Color',[0.90 0.45 0.10]);

xline(ax2, 0, '-', 'Color',[0.5 0.5 0.5]);
grid(ax2,'on');

leg2 = legend(ax2, [hGS hGD], {'Global Standard','Global Deviant'}, 'Location','southwest', FontSize=12);
set(leg2,'AutoUpdate','off');

hold(ax2,'off');
end

%% Figure mz4

function fig = plot_mz4_pupil(dataTableClean, timeVector, useField, participant, thisExp)
% plot_mz4_pupil
% 4-condition plot (Mean ± SEM):
%   Frequent Local Standard
%   Rare     Local Standard
%   Frequent Local Deviant
%   Rare     Local Deviant
%
% Usage (backward compatible):
%   fig = plot_mz4_pupil(dataTableClean, timeVector, 'PUPIL_Z')
%   fig = plot_mz4_pupil(dataTableClean, timeVector, 'PUPIL_Z', 'S01', 'E1')

% ---- Optional title prefix ----
titlePrefix = "";
if nargin >= 5 && ~isempty(participant) && ~isempty(thisExp)
    titlePrefix = sprintf("%s | %s | ", string(participant), string(thisExp));
end

fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w');
ax  = axes(fig); hold(ax,'on');
set(ax,'FontSize',12,'LineWidth',1);
set(ax,'Color','w','XColor','k','YColor','k');  % force white axes background
set(groot,'defaultAxesFontSize',11);  % affects all future axes
axis(ax,'tight');

title(ax, titlePrefix + "Pupil response by type", 'FontSize', 12);   % (Mean \pm SEM)
xlabel(ax, "Time (ms)", 'FontSize', 12);
ylabel(ax, "Pupil size (z-score)", 'FontSize', 12);

typeStr = string(dataTableClean.type);

condNames = [ ...
    "Frequent Local Standard", ...
    "Rare Local Standard", ...
    "Frequent Local Deviant", ...
    "Rare Local Deviant"];

% colors = [ ...
%     0.0000 0.4470 0.7410;
%     0.0000 0.4470 0.7410;
%     0.8500 0.3250 0.0980;
%     0.8500 0.3250 0.0980];

colors = [ ...
    0 0.45 0.85;
    0 0.45 0.85;
    0.90 0.45 0.10;
    0.90 0.45 0.10];

lineStyles = {'--','-','--','-'};

lineH = gobjects(numel(condNames), 1);

for k = 1:numel(condNames)
    thisCond = condNames(k);
    idxCond  = (typeStr == thisCond);
    if ~any(idxCond), continue; end

    pupMat = cell2mat(dataTableClean.(useField)(idxCond));
    meanCurve = mean(pupMat, 1, "omitnan");
    semCurve  = std(pupMat, 0, 1, "omitnan") ./ sqrt(size(pupMat,1));

    upper = meanCurve + semCurve;
    lower = meanCurve - semCurve;

    hFill = fill(ax, [timeVector fliplr(timeVector)], [upper fliplr(lower)], ...
                 colors(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    set(hFill, 'HandleVisibility','off');  % IMPORTANT: do not add fills to legend

    lineH(k) = plot(ax, timeVector, meanCurve, ...
        'Color', colors(k,:), ...
        'LineStyle', lineStyles{k}, ...
        'LineWidth', 2, ...
        'DisplayName', char(thisCond));
end

valid = isgraphics(lineH);
legend(ax, lineH(valid), 'Location', 'southwest', 'FontSize', 12);

grid(ax,'on'); box(ax,'off'); hold(ax,'off');

end

%% Figure Bekinschtein

function fig = plot_bekinschtein_pupil(dataTableClean, timeVector, useField, participant, thisExp)
% plot_bekinschtein_pupil
% Local Standard vs Local Deviant (Mean ± SEM)
%
% Usage (backward compatible):
%   fig = plot_bekinschtein_pupil(dataTableClean, timeVector, 'PUPIL_Z')
%   fig = plot_bekinschtein_pupil(dataTableClean, timeVector, 'PUPIL_Z', 'S01', 'E1')

% ---- Optional title prefix ----
titlePrefix = "";
if nargin >= 5 && ~isempty(participant) && ~isempty(thisExp)
    titlePrefix = sprintf("%s | %s | ", string(participant), string(thisExp));
end

fig = figure('Units','normalized','OuterPosition',[0 0 1 1], 'Color','w');
ax  = axes(fig); hold(ax,'on');
set(ax,'FontSize',12,'LineWidth',1);
set(ax,'Color','w','XColor','k','YColor','k');  % force white axes background
set(groot,'defaultAxesFontSize',11);  % affects all future axes
axis(ax,'tight');

title(ax, titlePrefix + "Pupil response by sound", 'FontSize', 12);   %(Mean \pm SEM)
xlabel(ax, "Time (ms)", 'FontSize', 12);
ylabel(ax, "Pupil size (z-score)", 'FontSize', 12);

soundStr   = string(dataTableClean.sound);
soundOrder = ["Local Standard","Local Deviant"];

% colors = [0 0.4470 0.7410;
%           0.8500 0.3250 0.0980];

colors = [0 0.45 0.85;
    0.90 0.45 0.10];

lineH = gobjects(numel(soundOrder),1);

for s = 1:numel(soundOrder)
    thisSound = soundOrder(s);
    idxSound  = (soundStr == thisSound);
    if ~any(idxSound), continue; end

    pupMat = cell2mat(dataTableClean.(useField)(idxSound));
    meanCurve = mean(pupMat, 1, "omitnan");
    semCurve  = std(pupMat, 0, 1, "omitnan") ./ sqrt(size(pupMat,1));

    upper = meanCurve + semCurve;
    lower = meanCurve - semCurve;

    hFill = fill(ax, [timeVector fliplr(timeVector)], [upper fliplr(lower)], ...
                 colors(s,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    set(hFill, 'HandleVisibility','off'); % IMPORTANT: do not add fills to legend

    lineH(s) = plot(ax, timeVector, meanCurve, ...
        'LineWidth', 2, 'Color', colors(s,:), ...
        'DisplayName', char(thisSound));
end

valid = isgraphics(lineH);
legend(ax, lineH(valid), 'Location', 'southwest', 'FontSize', 12);

grid(ax,'on'); box(ax,'off'); hold(ax,'off');

end


function d = prefer_exact_stem(d, exactName)
if isempty(d), return; end
idx=find(strcmpi({d.name},exactName),1);
if ~isempty(idx), d=d([idx setdiff(1:numel(d),idx,'stable')]); end
end

function fig=plot_ERtaVNS_visit_fourtrace(T,t,useField,participant,thisExp,visitNo,control)
% All retained phases, including habituation, as requested. Buffers NEVER
% enter either task averages or the task-reference normalization for control.
T=without_buffers(T);
visitLabel=get_visit_condition(T);
stimName="Sham"; if visitLabel=="ACTIVE", stimName="taVNS"; end
stim=valid_binary(T.stim_paired); pure=~stim;
local=string(T.sound); globalStatus=string(T.effect);
% Same standard/deviant colors; only line style distinguishes paired trials.
colors=[0 .447 .741; .85 .325 .098; 0 .447 .741; .85 .325 .098];
controlMatrix=[]; controlLabel="Pure "+stimName;
controlOverlay=struct('available',false,'source',control.source,'reason',control.reason);
if control.available
    assert(control.condition==visitLabel,'Control and task ACTIVE/SHAM metadata disagree.');
    assert(isequal(double(control.relTimes(:)'),double(t(:)')),'Control and task time axes differ.');
    % Existing PUPIL_Z is computed from all retained baseline-corrected task
    % samples. Reuse that exact reference for control in THIS figure; no
    % control-driven renormalization or changes to existing task curves.
    taskRaw=trial_matrix(T.PUPIL,numel(t));
    referenceMean=mean(taskRaw(:),'omitnan'); referenceSD=std(taskRaw(:),'omitnan');
    assert(isfinite(referenceSD) && referenceSD>0,'Task z-score reference is undefined.');
    taskZ=trial_matrix(T.(useField),numel(t)); expectedZ=(taskRaw-referenceMean)/referenceSD;
    delta=abs(taskZ-expectedZ); good=isfinite(delta);
    assert(any(good(:)) && max(delta(good))<1e-6,'Task normalization differs from current preprocess_pupil_LG; control overlay not rescaled by guessing.');
    controlMatrix=(trial_matrix(control.dataTableClean.PUPIL,numel(t))-referenceMean)/referenceSD;
    controlOverlay=struct('available',true,'source',control.source,'label',controlLabel, ...
        'n_retained',size(controlMatrix,1),'referenceMean',referenceMean,'referenceSD',referenceSD, ...
        'normalization',"task_reference_zscore",'anchor',"VNS_CONTROL intended stimulation onset");
end
fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 .8],'InvertHardcopy','off');
setappdata(fig,'controlOverlay',controlOverlay);
for panel=1:2
    ax=subplot(1,2,panel); hold(ax,'on');
    if panel==1
        standard=local=="Local Standard"; deviant=local=="Local Deviant";
        names=["Pure LS","Pure LD",stimName+" + LS",stimName+" + LD"]; panelName='Local';
    else
        standard=globalStatus=="Global Standard"; deviant=globalStatus=="Global Deviant";
        names=["Pure GS","Pure GD",stimName+" + GS",stimName+" + GD"]; panelName='Global';
    end
    masks=[pure & standard,pure & deviant,stim & standard,stim & deviant];
    styles={'-','-','--','--'}; h=gobjects(0);
    for k=1:4
        M=trial_matrix(T.(useField)(masks(:,k)),numel(t));
        if isempty(M), continue; end
        h(end+1)=plot_mean_sem(ax,t,M,names(k),styles{k},colors(k,:),2); %#ok<AGROW>
    end
    if ~isempty(controlMatrix)
        h(end+1)=plot_mean_sem(ax,t,controlMatrix,controlLabel,'-',[0 0 0],3.2); %#ok<AGROW>
    end
    xline(ax,0,'-','Color',[.5 .5 .5],'HandleVisibility','off');
    xlim(ax,[t(1) t(end)]); grid(ax,'on'); box(ax,'off');
    xlabel(ax,'Time from fifth tone / control stimulation onset (ms)'); ylabel(ax,'Pupil size (z-score)');
    title(ax,sprintf('%s | Visit %d | %s | %s | %s',participant,visitNo,thisExp,panelName,visitLabel));
    if ~isempty(h), legend(ax,h,'Location','best','Box','off','Interpreter','none'); end
end
sgtitle(fig,sprintf('%s | Visit %d | %s | Pure tones + %s paired trials + control | Mean +/- SEM', ...
    participant,visitNo,thisExp,stimName));
if ~control.available
    warning('C1:NoControlOverlay','%s V%d %s: only four task curves drawn; %s.',participant,visitNo,thisExp,control.reason);
    sgtitle(fig,sprintf('%s | Visit %d | %s | Pure tones + %s paired trials | Mean +/- SEM | Control unavailable', ...
        participant,visitNo,thisExp,stimName));
end
end

function h=plot_mean_sem(ax,t,M,label,style,color,width)
n=size(M,1); mu=mean(M,1,'omitnan'); finiteN=sum(isfinite(M),1);
sem=std(M,0,1,'omitnan')./sqrt(finiteN); sem(finiteN<2)=NaN;
fill(ax,[t fliplr(t)],[mu+sem fliplr(mu-sem)],color, ...
    'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off','Clipping','on');
h=plot(ax,t,mu,'LineStyle',style,'Color',color,'LineWidth',width, ...
    'DisplayName',sprintf('%s (n=%d)',label,n),'Clipping','on');
end

function apply_pupil_limits(fig,upperLimit)
axesList=findall(fig,'Type','axes'); lower=Inf;
for k=1:numel(axesList), limits=ylim(axesList(k)); lower=min(lower,limits(1)); end
if ~isfinite(lower) || lower>=upperLimit, lower=upperLimit-1; end
for k=1:numel(axesList)
    ax=axesList(k); ylim(ax,[lower upperLimit]); set(ax,'Clipping','on');
    yticks(ax,unique([yticks(ax) upperLimit]));
end
end

function M=trial_matrix(C,nt)
if isempty(C), M=zeros(0,nt); return; end
assert(iscell(C) && all(cellfun(@numel,C)==nt),'Unexpected pupil trial matrix.');
M=nan(numel(C),nt);
for k=1:numel(C), M(k,:)=double(C{k}(:)'); end
end

function T=without_buffers(T)
assert(ismember('is_post_stim_buffer',T.Properties.VariableNames),'Missing buffer metadata; rerun current preprocessing.');
T=T(~valid_binary(T.is_post_stim_buffer),:);
end

function b=valid_binary(x)
assert(isnumeric(x) || islogical(x),'Boolean metadata must be numeric/logical.');
assert(all(isfinite(x(:)) & (x(:)==0 | x(:)==1)),'Missing/invalid boolean metadata.');
b=logical(x(:));
end

function condition=get_visit_condition(T)
assert(ismember('session_condition',T.Properties.VariableNames),'Missing ACTIVE/SHAM metadata.');
c=upper(strtrim(string(T.session_condition)));
c=c(~ismissing(c) & c~=""); c=unique(c);
assert(isscalar(c) && ismember(c,["ACTIVE","SHAM"]),'Missing or conflicting ACTIVE/SHAM metadata.');
condition=c;
end

function lg=verify_buffer_sequence(lg)
T=lg.trialInfo; needed={'trial','block','global_short','stim_paired','t5_time_ms','is_post_stim_buffer'};
assert(all(ismember(needed,T.Properties.VariableNames)),'Incomplete raw trial metadata; use current parse_LG_eyelink.');
stim=valid_binary(T.stim_paired); candidates=false(height(T),1); violations=zeros(0,1);
assert(isequal(double(T.trial(:)),(1:height(T))'),'Raw trial identifiers are not the expected indexed sequence.');
for k=find(stim)'
    next=k+1;
    valid=next<=height(T);
    if valid
        valid=isfinite(T.t5_time_ms(k)) && isfinite(T.t5_time_ms(next)) && ...
            isfinite(T.block(k)) && T.block(next)==T.block(k) && ...
            upper(string(T.global_short(next)))=="GS" && ~stim(next);
    end
    if valid, candidates(next)=true; else, violations(end+1,1)=k; end %#ok<AGROW>
end
detected=any(stim) && isempty(violations);
% Protect pre-buffer recordings: never exclude accidental N+1 trials when
% the full mandatory rule does not hold. Existing parser uses the same rule.
flags=candidates & detected;
lg.trialInfo.is_post_stim_buffer=flags;
if ismember('exclude_primary_pupil',T.Properties.VariableNames), lg.trialInfo.exclude_primary_pupil=flags; end
if ismember('exclude_primary_eeg',T.Properties.VariableNames), lg.trialInfo.exclude_primary_eeg=flags; end
lg.report.bufferRuleDetected=detected; lg.report.nPostStimBuffers=nnz(flags);
lg.report.nBufferRuleViolations=numel(violations); lg.report.bufferViolationTrials=violations;
fprintf('[C1 buffer check] rule verified=%d | buffers=%d | violations=%d\n',detected,nnz(flags),numel(violations));
if any(stim) && ~detected
    warning('C1:BufferRuleNotMet','Complete N+1 rule not met. No accidental following GS trials are reclassified; inspect bufferQC.');
end
end

function assert_buffer_exclusion(T,raw)
assert(ismember('is_post_stim_buffer',T.Properties.VariableNames),'Missing clean buffer flags.');
assert(~any(valid_binary(T.is_post_stim_buffer)),'A buffer survived preprocessing; plots stopped.');
assert(ismember('trialInFile',T.Properties.VariableNames),'Cannot cross-check clean trial identities.');
idx=double(T.trialInFile);
assert(all(isfinite(idx) & idx>=1 & idx<=height(raw) & idx==fix(idx)),'Invalid clean trial identities.');
assert(~any(raw.is_post_stim_buffer(idx)),'A raw-sequence buffer survived preprocessing; plots stopped.');
end

function C=preprocess_control_pupil(lines,visitNo,source,fs,dt,pad_ms,tout,b_lp,a_lp,blinkLimit,baseWin)
% Standalone stimulation-only processing; NO invented tone or GS/GD labels.
% Same continuous blink interpolation/filter, epoch interpolation, baseline,
% blink threshold and incomplete-epoch rejection as preprocess_pupil_LG.
assert(numel(baseWin)==2 && all(isfinite(baseWin)) && baseWin(1)<baseWin(2) && baseWin(2)<=0,'Invalid control baseline.');
assert(abs(dt-1000/fs)<1e-9,'Pupil sampling settings disagree.');
% A0 uses ERtaVNS_condition_site: "left cymba conchae" / "left earlobe".
% SITE therefore contains spaces. End that field at VISIT_, not at a space.
pat=['^MSG\s+([0-9]+(?:\.[0-9]+)?)\s+VNS_CONTROL\s+BURST_(\d+)\s+' ...
    'CODE_(\d+)\s+CONDITION_(ACTIVE|SHAM)\s+SITE_(.+?)\s+VISIT_(\d+)\s*$'];
onsets=[]; bursts=[]; conditions=strings(0,1); sites=strings(0,1); visits=[];
sampleTimes=[]; samplePupil=[]; blinkStart=[]; blinkEnd=[];
for k=1:numel(lines)
    L=strtrim(lines{k}); if isempty(L), continue; end
    tok=regexp(L,pat,'tokens','once');
    if ~isempty(tok)
        assert(str2double(tok{3})==10,'Unexpected control marker code.');
        onsets(end+1,1)=str2double(tok{1}); bursts(end+1,1)=str2double(tok{2}); %#ok<AGROW>
        conditions(end+1,1)=string(tok{4}); visits(end+1,1)=str2double(tok{6}); %#ok<AGROW>
        sites(end+1,1)=strtrim(string(tok{5})); %#ok<AGROW>
    elseif ~isempty(regexp(L,'^MSG\s+\S+\s+VNS_CONTROL(?:\s|$)','once'))
        error('C1:ControlMarkerFormat','Unrecognized control onset message: %s',L);
    end
    if isstrprop(L(1),'digit')
        vals=sscanf(L,'%f',5);
        if numel(vals)>=5
            sampleTimes(end+1,1)=vals(1); samplePupil(end+1,1)=vals(4); %#ok<AGROW>
        end
    end
    tok=regexp(L,'^EBLINK\s+\w+\s+(\d+)\s+(\d+)\s+(\d+)','tokens','once');
    if ~isempty(tok)
        blinkStart(end+1,1)=str2double(tok{1}); blinkEnd(end+1,1)=str2double(tok{2}); %#ok<AGROW>
    end
end
assert(~isempty(onsets),'No VNS_CONTROL timestamps found. Convert the same-visit control EDF to ASC.');
assert(all(visits==visitNo),'Control marker visit disagrees with folder visit.');
condition=unique(conditions); assert(isscalar(condition),'Conflicting control ACTIVE/SHAM metadata.');
assert(numel(unique(bursts))==numel(bursts) && all(diff(onsets)>0),'Duplicate/nonchronological control markers.');
assert(~isempty(sampleTimes),'No pupil samples parsed in control ASC.');
if numel(onsets)~=30, warning('C1:ControlCount','Control contains %d burst markers (design plans 30); using observed markers only.',numel(onsets)); end
blink=false(size(sampleTimes));
for k=1:numel(blinkStart)
    blink(sampleTimes>=blinkStart(k)-pad_ms & sampleTimes<=blinkEnd(k)+pad_ms)=true;
end
pupil=samplePupil; pupil(blink)=NaN;
pupil=fillmissing(pupil,'linear'); pupil=filtfilt(b_lp,a_lp,pupil);
pre=dt*floor((min(tout(1),baseWin(1))-100)/dt);
post=dt*ceil((max(tout(end),baseWin(2))+200)/dt); big=pre:dt:post;
outIdx=round((tout-pre)/dt)+1; baseIdx=big>=baseWin(1) & big<=baseWin(2);
assert(nnz(baseIdx)>=5,'Fewer than five control baseline samples.');
N=numel(onsets); epochs=cell(N,1); blinkEpoch=cell(N,1); reason=strings(N,1);
for k=1:N
    if onsets(k)+pre<sampleTimes(1) || onsets(k)+post>sampleTimes(end)
        reason(k)="Not fully inside recording window"; continue;
    end
    idx=sampleTimes>=onsets(k)+pre & sampleTimes<=onsets(k)+post;
    ts=sampleTimes(idx)-onsets(k); pp=pupil(idx); bb=blink(idx);
    good=isfinite(ts); ts=ts(good); pp=pp(good); bb=bb(good);
    [ts,order]=sort(ts); pp=pp(order); bb=bb(order);
    [ts,uniqueIdx]=unique(ts,'stable'); pp=pp(uniqueIdx); bb=bb(uniqueIdx);
    if numel(ts)<2, reason(k)="No samples in epoch window"; continue; end
    y=interp1(ts,pp,big,'linear',NaN);
    b=interp1(ts,double(bb),big,'nearest',0)>0;
    y=y-mean(y(baseIdx),'omitnan'); y=y(outIdx); b=b(outIdx);
    if mean(b)>blinkLimit, reason(k)="Blink threshold exceeded"; continue; end
    if any(~isfinite(y)), reason(k)="Trial incomplete (missing data)"; continue; end
    epochs{k}=y; blinkEpoch{k}=b;
end
keep=reason=="";
T=table(bursts(keep),onsets(keep),conditions(keep),visits(keep),sites(keep),epochs(keep),blinkEpoch(keep), ...
    'VariableNames',{'burst','time','session_condition','visit','stim_site','PUPIL','BLINK'});
M=trial_matrix(T.PUPIL,numel(tout)); zMean=NaN; zSD=NaN;
if ~isempty(M)
    zMean=mean(M(:),'omitnan'); zSD=std(M(:),'omitnan');
    if isfinite(zSD) && zSD>0, Z=(M-zMean)/zSD; else, Z=nan(size(M)); end
    T.PUPIL_Z=num2cell(Z,2);
else
    T.PUPIL_Z=cell(0,1);
end
removed=table(bursts(~keep),onsets(~keep),reason(~keep),'VariableNames',{'burst','time','Reason'});
qc=struct('source',string(source),'visit',visitNo,'condition',condition,'n_bursts',N,'n_retained',nnz(keep), ...
    'n_rejected',nnz(~keep),'anchor',"VNS_CONTROL intended stimulation onset", ...
    'baseline_ms',baseWin,'epoch_ms',[tout(1) tout(end)],'fs',fs,'blink_padding_ms',pad_ms, ...
    'blink_limit',blinkLimit,'normalization_mean',zMean,'normalization_sd',zSD, ...
    'overlay_normalization',"Uses the matched task reference; stored PUPIL_Z is control-only normalization.");
why=""; if ~any(keep), why="No control bursts survived preprocessing"; end
C=struct('available',any(keep),'condition',condition,'dataTableClean',T,'source',string(source), ...
    'reason',why,'removedTrials',removed,'qc',qc,'relTimes',tout);
end
