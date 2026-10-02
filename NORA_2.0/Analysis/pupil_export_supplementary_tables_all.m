function OUT = pupil_export_supplementary_tables_all(outpath, grandOutFolder, STATS, GRAND, SSI, Tinteraction, OUT_PUPIL_RAINBOX, PUPIL_GROUP_DESC, excludeSubs)

tablesDir = fullfile(grandOutFolder, 'tables');
if ~exist(tablesDir, 'dir')
    mkdir(tablesDir);
end

if nargin < 9 || isempty(excludeSubs)
    excludeSubs = strings(0,1);
end
excludeSubs = string(excludeSubs);

expList   = ["E1","E2"];
expLabels = ["Counting","Visual interference"];

% ============================================================
% Load preprocessing summary table (the one with:
% Num_bad_trials_E1_Pupil, Num_bad_trials_E2_Pupil,
% Num_good_trials_E1_Pupil, Num_good_trials_E2_Pupil)
% ============================================================
PQA = pupil_try_load_info_table(outpath, grandOutFolder);

%% ============================================================
% S1a / S1b: subject-level preprocessing / retention
% ============================================================
TsampleLong = table();

d = dir(fullfile(outpath, 'S*'));
d = d([d.isdir]);
subNames = string({d.name});
subNames = subNames(~ismember(subNames, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(subNames), '^S\d+$', 'once'));
subNames = subNames(keep);
subNames = subNames(~ismember(subNames, excludeSubs));

for i = 1:numel(expList)
    E = expList(i);

    for s = 1:numel(subNames)
        participant = subNames(s);

        matFile = fullfile(outpath, char(participant), sprintf('%s_%s_pupil_preproc.mat', participant, E));
        if ~exist(matFile, 'file')
            continue;
        end

        S = load(matFile, 'dataTableClean');
        if ~isfield(S,'dataTableClean') || isempty(S.dataTableClean) || height(S.dataTableClean)==0
            continue;
        end

        T = S.dataTableClean;

        % retained trial counts from cleaned pupil table
        soundVec = string(T.sound);
        typeVec  = string(T.type);

        nLS  = sum(soundVec == "Local Standard");
        nLD  = sum(soundVec == "Local Deviant");
        nFLS = sum(typeVec == "Frequent Local Standard");
        nRLS = sum(typeVec == "Rare Local Standard");
        nFLD = sum(typeVec == "Frequent Local Deviant");
        nRLD = sum(typeVec == "Rare Local Deviant");

        effectVec = strings(height(T),1);
        effectVec(startsWith(typeVec, "Frequent")) = "Global Standard";
        effectVec(startsWith(typeVec, "Rare"))     = "Global Deviant";

        nGS = sum(effectVec == "Global Standard");
        nGD = sum(effectVec == "Global Deviant");

        % preprocessing totals from info table
        qa = pupil_lookup_preproc_row(PQA, participant, E);

        TsampleLong = [TsampleLong; table( ...
            participant, expLabels(i), E, ...
            qa.GoodTrials, qa.RejectedTrials, ...
            nLS, nLD, nFLS, nRLS, nFLD, nRLD, nGS, nGD, ...
            'VariableNames', { ...
            'Subject','ExperimentLabel','ExperimentKey','GoodTrials_Total','RejectedTrials_Total', ...
            'Retained_LS','Retained_LD','Retained_FLS','Retained_RLS','Retained_FLD','Retained_RLD','Retained_GS','Retained_GD'})];
    end
end

Tsample = table();
for i = 1:numel(expList)
    idx = TsampleLong.ExperimentKey == expList(i);
    if ~any(idx), continue; end

    [m0,se0] = pupil_mean_se(TsampleLong.GoodTrials_Total(idx));
    [mr,ser] = pupil_mean_se(TsampleLong.RejectedTrials_Total(idx));
    [m1,se1] = pupil_mean_se(TsampleLong.Retained_LS(idx));
    [m2,se2] = pupil_mean_se(TsampleLong.Retained_LD(idx));
    [m3,se3] = pupil_mean_se(TsampleLong.Retained_FLS(idx));
    [m4,se4] = pupil_mean_se(TsampleLong.Retained_RLS(idx));
    [m5,se5] = pupil_mean_se(TsampleLong.Retained_FLD(idx));
    [m6,se6] = pupil_mean_se(TsampleLong.Retained_RLD(idx));
    [m7,se7] = pupil_mean_se(TsampleLong.Retained_GS(idx));
    [m8,se8] = pupil_mean_se(TsampleLong.Retained_GD(idx));

    Tsample = [Tsample; table( ...
        expLabels(i), expList(i), sum(idx), ...
        m0,se0, mr,ser, ...
        m1,se1, m2,se2, ...
        m3,se3, m4,se4, m5,se5, m6,se6, ...
        m7,se7, m8,se8, ...
        'VariableNames', { ...
        'Experiment','ExperimentKey','Nsub', ...
        'MeanGoodTrials_Total','SEGoodTrials_Total', ...
        'MeanRejectedTrials_Total','SERejectedTrials_Total', ...
        'MeanRetained_LS','SERetained_LS', ...
        'MeanRetained_LD','SERetained_LD', ...
        'MeanRetained_FLS','SERetained_FLS', ...
        'MeanRetained_RLS','SERetained_RLS', ...
        'MeanRetained_FLD','SERetained_FLD', ...
        'MeanRetained_RLD','SERetained_RLD', ...
        'MeanRetained_GS','SERetained_GS', ...
        'MeanRetained_GD','SERetained_GD'})];
end

%% ============================================================
% S2: Group-level cluster summary
% ============================================================
Tgroup = table();

for i = 1:numel(expList)
    E = expList(i);
    t = STATS.(E).relTimes(:)';

    % Local
    mask = STATS.(E).maskLocal(:)';
    dvec = GRAND.(E).mu_sound(2,:) - GRAND.(E).mu_sound(1,:);
    cls  = STATS.(E).clustersLocal;

    if any(mask)
        dd = diff([false mask false]);
        rs = find(dd==1);
        re = find(dd==-1)-1;

        pBest = min([cls.pFWE]);
        idx = find(mask);
        [~,ix] = max(abs(dvec(idx)));
        idxPeak = idx(ix);

        row = table( ...
            expLabels(i), "Local", ...
            t(rs(1)), t(re(end)), ...
            pBest, dvec(idxPeak), t(idxPeak), ...
            string(sign(dvec(idxPeak))>0), ...
            'VariableNames', {'Experiment','Effect','StartMS','EndMS','Best_pFWE','PeakDiff','PeakMS','PositiveDirection'});
    else
        row = table( ...
            expLabels(i), "Local", ...
            NaN, NaN, NaN, NaN, NaN, "false", ...
            'VariableNames', {'Experiment','Effect','StartMS','EndMS','Best_pFWE','PeakDiff','PeakMS','PositiveDirection'});
    end
    Tgroup = [Tgroup; row];

    % Global
    mask = STATS.(E).maskGlobal(:)';
    dvec = GRAND.(E).mu_effect(2,:) - GRAND.(E).mu_effect(1,:);
    cls  = STATS.(E).clustersGlobal;

    if any(mask)
        dd = diff([false mask false]);
        rs = find(dd==1);
        re = find(dd==-1)-1;

        pBest = min([cls.pFWE]);
        idx = find(mask);
        [~,ix] = max(abs(dvec(idx)));
        idxPeak = idx(ix);

        row = table( ...
            expLabels(i), "Global", ...
            t(rs(1)), t(re(end)), ...
            pBest, dvec(idxPeak), t(idxPeak), ...
            string(sign(dvec(idxPeak))>0), ...
            'VariableNames', {'Experiment','Effect','StartMS','EndMS','Best_pFWE','PeakDiff','PeakMS','PositiveDirection'});
    else
        row = table( ...
            expLabels(i), "Global", ...
            NaN, NaN, NaN, NaN, NaN, "false", ...
            'VariableNames', {'Experiment','Effect','StartMS','EndMS','Best_pFWE','PeakDiff','PeakMS','PositiveDirection'});
    end
    Tgroup = [Tgroup; row];
end

%% ============================================================
% S2b: Group-level cluster descriptives
% ============================================================
if nargin < 8 || isempty(PUPIL_GROUP_DESC)
    TgroupDesc = table();
else
    TgroupDesc = PUPIL_GROUP_DESC.table;
end

%% ============================================================
% S3 / S3b: Secondary between-task comparison
% ============================================================
Train = OUT_PUPIL_RAINBOX.statsTable;
Tsub  = OUT_PUPIL_RAINBOX.subjectTable;

%% ============================================================
% S4: Individual-level SSI summary
% ============================================================
Tind = table();

for i = 1:numel(expList)
    E = expList(i);
    t = SSI.(E).times(:)';
    keep = t >= 0;
    t2 = t(keep);
    dt = median(diff(t2));

    for eff = ["Local","Global"]
        if eff == "Local"
            M = SSI.(E).local(:,keep);
        else
            M = SSI.(E).global(:,keep);
        end

        hasEff = any(M > 0, 2);
        dur = sum(M > 0, 2) * dt;
        meanSSI = mean(M,1,'omitnan');
        [peakSSI,ix] = max(meanSSI);

        row = table( ...
            expLabels(i), eff, ...
            sum(hasEff), size(M,1), 100*sum(hasEff)/size(M,1), ...
            peakSSI, t2(ix), ...
            median(dur,'omitnan'), ...
            median(dur(hasEff),'omitnan'), ...
            'VariableNames', {'Experiment','Effect','N_with_effect','N_total','Percent_with_effect', ...
                              'PeakMeanSSI','PeakMS','MedianDuration_All_ms','MedianDuration_EffectOnly_ms'});
        Tind = [Tind; row];
    end
end

%% ============================================================
% S5: Interaction summary
% ============================================================
Tint = Tinteraction;

%% ============================================================
% Save CSV
% ============================================================
writetable(TsampleLong, fullfile(tablesDir, 'PUPIL_S1a_Preprocessing_Long.csv'));
writetable(Tsample,     fullfile(tablesDir, 'PUPIL_S1b_Preprocessing_Summary.csv'));
writetable(Tgroup,      fullfile(tablesDir, 'PUPIL_S2_Primary_GroupEffects.csv'));
writetable(TgroupDesc,  fullfile(tablesDir, 'PUPIL_S2b_Primary_GroupDescriptives.csv'));
writetable(Train,       fullfile(tablesDir, 'PUPIL_S3_SecondaryComparisons.csv'));
writetable(Tsub,        fullfile(tablesDir, 'PUPIL_S3b_Secondary_SubjectValues.csv'));
writetable(Tind,        fullfile(tablesDir, 'PUPIL_S4_IndividualSSI.csv'));
writetable(Tint,        fullfile(tablesDir, 'PUPIL_S5_InteractionSummary.csv'));

%% ============================================================
% Save Excel
% ============================================================
xlsxFile = fullfile(tablesDir, 'PUPIL_Supplementary_Tables.xlsx');
writetable(TsampleLong, xlsxFile, 'Sheet', 'S1a_Preproc_Long');
writetable(Tsample,     xlsxFile, 'Sheet', 'S1b_Preproc_Sum');
writetable(Tgroup,      xlsxFile, 'Sheet', 'S2_Group');
writetable(TgroupDesc,  xlsxFile, 'Sheet', 'S2b_GroupDesc');
writetable(Train,       xlsxFile, 'Sheet', 'S3_Secondary');
writetable(Tsub,        xlsxFile, 'Sheet', 'S3b_SubjectVals');
writetable(Tind,        xlsxFile, 'Sheet', 'S4_IndividualSSI');
writetable(Tint,        xlsxFile, 'Sheet', 'S5_Interaction');

OUT = struct();
OUT.tablesDir = tablesDir;
OUT.sampleLongTable = TsampleLong;
OUT.sampleTable = Tsample;
OUT.groupTable = Tgroup;
OUT.groupDescTable = TgroupDesc;
OUT.rainboxStatsTable = Train;
OUT.rainboxSubjectTable = Tsub;
OUT.individualTable = Tind;
OUT.interactionTable = Tint;

fprintf('\n[PUPIL TABLES] Saved supplementary tables to:\n%s\n', tablesDir);
end

% ============================================================
function [mu,se] = pupil_mean_se(x)
x = double(x);
x = x(isfinite(x));
if isempty(x)
    mu = NaN;
    se = NaN;
else
    mu = mean(x,'omitnan');
    se = std(x,0,'omitnan') / sqrt(numel(x));
end
end

% ============================================================
function PQA = pupil_try_load_info_table(outpath, grandOutFolder)

PQA = struct();
PQA.exists = false;
PQA.T = table();
PQA.source = "";

candidates = {
    fullfile(outpath, 'Info_EEG_LG.mat')
    fullfile(grandOutFolder, 'Info_EEG_LG.mat')
    fullfile(fileparts(grandOutFolder), 'Info_EEG_LG.mat')
};

for i = 1:numel(candidates)
    f = candidates{i};
    if exist(f, 'file') ~= 2
        continue;
    end

    S = load(f);
    fn = string(fieldnames(S));

    if ismember("T_LG", fn) && istable(S.T_LG)
        PQA.exists = true;
        PQA.T = S.T_LG;
        PQA.source = f;
        return;
    end

    for k = 1:numel(fn)
        if istable(S.(fn(k)))
            PQA.exists = true;
            PQA.T = S.(fn(k));
            PQA.source = f;
            return;
        end
    end
end
end

% ============================================================
function row = pupil_lookup_preproc_row(PQA, subjectID, expKey)

row = struct();
row.GoodTrials = NaN;
row.RejectedTrials = NaN;

if ~isfield(PQA,'exists') || ~PQA.exists || isempty(PQA.T)
    return;
end

T = PQA.T;
vn = string(T.Properties.VariableNames);

% subject column
subCol = "";
for cand = ["Subject","subject","Sub","sub","Subject_ID","participant","Participant"]
    if ismember(cand, vn)
        subCol = cand;
        break;
    end
end

if subCol == ""
    return;
end

subVals = string(T.(subCol));
m = strcmpi(strtrim(subVals), string(subjectID));
if ~any(m)
    m = strcmpi(erase(strtrim(subVals)," "), erase(strtrim(string(subjectID))," "));
end
if ~any(m)
    return;
end

r = T(find(m,1),:);

if expKey == "E1"
    if ismember("Num_bad_trials_E1_Pupil", vn)
        row.RejectedTrials = double(r.("Num_bad_trials_E1_Pupil"));
    end
    if ismember("Num_good_trials_E1_Pupil", vn)
        row.GoodTrials = double(r.("Num_good_trials_E1_Pupil"));
    end
elseif expKey == "E2"
    if ismember("Num_bad_trials_E2_Pupil", vn)
        row.RejectedTrials = double(r.("Num_bad_trials_E2_Pupil"));
    end
    if ismember("Num_good_trials_E2_Pupil", vn)
        row.GoodTrials = double(r.("Num_good_trials_E2_Pupil"));
    end
end
end