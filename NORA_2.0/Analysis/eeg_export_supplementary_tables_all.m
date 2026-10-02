function OUT = eeg_export_supplementary_tables_all( ...
    outroot, grandFold, GROUPPMAP, SSI, cfg, excludeSubs, targets, roiChans, OUT_EEG_RAINBOX_PZ, INT)

%% SUMMARY TABLES INFORMATION
% ==============================================================
% eeg_export_supplementary_tables_all
%
% Exports EEG supplementary tables to:
%   GRAND_AVG\tables\
%
% Outputs:
%   S1  Sample / task / retention
%   S2  Group scalp-level effects
%   S3  ROI follow-up table (Pz, Fz)
%   S4  Individual-level SSI summary
%   S5  Secondary Pz comparison table (from rainbox output)
%
% Saves:
%   - one Excel workbook (.xlsx)
%   - one Word document (.docx)
%   - separate CSV files
%
% NOTE:
%   The S1 QA columns (rejected components / trials / channels) try to
%   auto-read Info_EEG_LG.mat if available. If your variable names differ,
%   those columns may appear as NaN and may need a tiny adjustment.
% ==============================================================

if nargin < 6 || isempty(excludeSubs)
    excludeSubs = strings(0,1);
end
if nargin < 7 || isempty(targets)
    targets = [ ...
        "Frequent Local Standard", ...
        "Rare Local Standard", ...
        "Frequent Local Deviant", ...
        "Rare Local Deviant"];
end
if nargin < 8 || isempty(roiChans)
    roiChans = ["Pz","Fz"];
end
if nargin < 9
    OUT_EEG_RAINBOX_PZ = [];
end

if nargin < 10
    INT = struct();
end

excludeSubs = string(excludeSubs);
targets     = string(targets);
roiChans    = string(roiChans);

if ~exist(grandFold, 'dir')
    mkdir(grandFold);
end

tablesDir = fullfile(grandFold, 'tables');
if ~exist(tablesDir, 'dir')
    mkdir(tablesDir);
end

expList   = ["E1","E2"];
expLabels = ["Counting","Visual interference"];
effList   = ["Local","Global"];

%% ------------------------------------------------------------
% Try to load QA / preprocessing info table
%% ------------------------------------------------------------
QA = supp2_load_info_table(outroot, grandFold);

%% ------------------------------------------------------------
% S1: SAMPLE / TASK / RETENTION
%% ------------------------------------------------------------
TretLong = table();

for ei = 1:numel(expList)
    E = expList(ei);

    subs = supp2_list_subjects_with_file(outroot, E, excludeSubs);
    for si = 1:numel(subs)
        sID = subs(si);

        setFile = fullfile(outroot, char(sID), sprintf('auditory_all_LG_%s_%s.set', char(sID), char(E)));
        if ~exist(setFile,'file')
            continue;
        end

        EEG = pop_loadset(setFile);
        EEG = eeg_checkset(EEG);

        counts = supp2_count_trials_from_labels(EEG, targets);
        qaRow   = supp2_lookup_qa_row(QA, sID, E);

        TretLong = [TretLong; table( ...
            sID, expLabels(ei), E, ...
            EEG.trials, ...
            counts.FLS, counts.RLS, counts.FLD, counts.RLD, ...
            counts.LocalStandard, counts.LocalDeviant, ...
            counts.GlobalStandard, counts.GlobalDeviant, ...
            qaRow.RejectedComponents, qaRow.RejectedTrials, ...
            qaRow.GoodTrials, qaRow.RejectedChannels, ...
            qaRow.InfoSource, ...
            'VariableNames', { ...
            'Subject','ExperimentLabel','ExperimentKey', ...
            'GoodTrials_ALLset', ...
            'Retained_FLS','Retained_RLS','Retained_FLD','Retained_RLD', ...
            'LocalStandardTrials','LocalDeviantTrials', ...
            'GlobalStandardTrials','GlobalDeviantTrials', ...
            'EEG_RejectedComponents','EEG_RejectedTrials','EEG_GoodTrials_QA','EEG_RejectedChannels', ...
            'InfoSource'})];
    end
end

TretSummary = table();
for ei = 1:numel(expList)
    idx = TretLong.ExperimentKey == expList(ei);
    if ~any(idx), continue; end

    [m1,se1] = supp2_mean_se(TretLong.GoodTrials_ALLset(idx));
    [m2,se2] = supp2_mean_se(TretLong.LocalStandardTrials(idx));
    [m3,se3] = supp2_mean_se(TretLong.LocalDeviantTrials(idx));
    [m4,se4] = supp2_mean_se(TretLong.GlobalStandardTrials(idx));
    [m5,se5] = supp2_mean_se(TretLong.GlobalDeviantTrials(idx));
    [m6,se6] = supp2_mean_se(TretLong.EEG_RejectedComponents(idx));
    [m7,se7] = supp2_mean_se(TretLong.EEG_RejectedTrials(idx));
    [m8,se8] = supp2_mean_se(TretLong.EEG_GoodTrials_QA(idx));
    [m9,se9] = supp2_mean_se(TretLong.EEG_RejectedChannels(idx));

    TretSummary = [TretSummary; table( ...
        expLabels(ei), expList(ei), ...
        sum(idx), ...
        m1,se1, ...
        m2,se2, ...
        m3,se3, ...
        m4,se4, ...
        m5,se5, ...
        m6,se6, ...
        m7,se7, ...
        m8,se8, ...
        m9,se9, ...
        'VariableNames', { ...
        'ExperimentLabel','ExperimentKey','Nsubjects', ...
        'MeanGoodTrials_ALLset','SEGoodTrials_ALLset', ...
        'MeanLocalStandardTrials','SELocalStandardTrials', ...
        'MeanLocalDeviantTrials','SELocalDeviantTrials', ...
        'MeanGlobalStandardTrials','SEGlobalStandardTrials', ...
        'MeanGlobalDeviantTrials','SEGlobalDeviantTrials', ...
        'MeanRejectedComponents','SERejectedComponents', ...
        'MeanRejectedTrials','SERejectedTrials', ...
        'MeanGoodTrials_QA','SEGoodTrials_QA', ...
        'MeanRejectedChannels','SERejectedChannels'})];
end

%% ------------------------------------------------------------
% S2: GROUP SCALP-LEVEL EFFECTS
%% ------------------------------------------------------------
Tgroup = table();

for ei = 1:numel(expList)
    E = expList(ei);
    if ~isfield(GROUPPMAP, E)
        continue;
    end

    chanLabels = supp2_get_template_chanlabels(outroot, E, excludeSubs);

    for fi = 1:numel(effList)
        effectName = effList(fi);

        if effectName == "Local"
            pCT = GROUPPMAP.(E).pLoc;
            m   = logical(GROUPPMAP.(E).mLoc(:)');
        else
            pCT = GROUPPMAP.(E).pGlo;
            m   = logical(GROUPPMAP.(E).mGlo(:)');
        end

        t = double(GROUPPMAP.(E).times(:)');
        [rs, re] = supp2_find_runs(m);

        if isempty(rs)
            Tgroup = [Tgroup; table( ...
                expLabels(ei), E, effectName, false, ...
                NaN, NaN, NaN, NaN, NaN, NaN, "", ...
                string(cfg.mcc.method), cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp, ...
                'VariableNames', { ...
                'ExperimentLabel','ExperimentKey','Effect','HasSignificantWindow', ...
                'WindowIndex','StartMS','EndMS','DurationMS','MinPFWE','MinPTimeMS','MinPChannel', ...
                'MCC_Method','Alpha_Group','MinElectrodes','MinConsecutiveSamples'})];
            continue;
        end

        for w = 1:numel(rs)
            idxWin = rs(w):re(w);
            pWin   = pCT(:, idxWin);

            [pmin, linIdx] = min(pWin(:));
            [chIdx, tLocalIdx] = ind2sub(size(pWin), linIdx);
            tMin = t(idxWin(tLocalIdx));

            if isempty(chanLabels) || chIdx > numel(chanLabels)
                chName = "";
            else
                chName = string(chanLabels{chIdx});
            end

            Tgroup = [Tgroup; table( ...
                expLabels(ei), E, effectName, true, ...
                w, t(rs(w)), t(re(w)), t(re(w)) - t(rs(w)), ...
                pmin, tMin, chName, ...
                string(cfg.mcc.method), cfg.alpha_grp, cfg.minElec_grp, cfg.minConsec_grp, ...
                'VariableNames', { ...
                'ExperimentLabel','ExperimentKey','Effect','HasSignificantWindow', ...
                'WindowIndex','StartMS','EndMS','DurationMS','MinPFWE','MinPTimeMS','MinPChannel', ...
                'MCC_Method','Alpha_Group','MinElectrodes','MinConsecutiveSamples'})];
        end
    end
end

%% ------------------------------------------------------------
% S3: ROI FOLLOW-UP TABLE (Pz / Fz)
%% ------------------------------------------------------------
Troi = table();

for ei = 1:numel(expList)
    E = expList(ei);
    if ~isfield(GROUPPMAP, E)
        continue;
    end

    subs = supp2_list_subjects_with_file(outroot, E, excludeSubs);
    if isempty(subs)
        continue;
    end

    chanLabels = supp2_get_template_chanlabels(outroot, E, excludeSubs);

    roiData = struct();
    for ri = 1:numel(roiChans)
        roi = roiChans(ri);
        roiField = matlab.lang.makeValidName(char(roi));

        roiData.(roiField).subject = strings(0,1);
        roiData.(roiField).C = {};

        for si = 1:numel(subs)
            sID = subs(si);
            setFile = fullfile(outroot, char(sID), sprintf('auditory_all_LG_%s_%s.set', char(sID), char(E)));
            if ~exist(setFile,'file')
                continue;
            end

            C = supp2_load_subject_curves_eeg(setFile, roi, targets);
            if isempty(C)
                continue;
            end

            roiData.(roiField).subject(end+1,1) = sID; %#ok<AGROW>
            roiData.(roiField).C{end+1,1} = C; %#ok<AGROW>
        end
    end

    for fi = 1:numel(effList)
        effectName = effList(fi);

        if effectName == "Local"
            pCT = GROUPPMAP.(E).pLoc;
            m   = logical(GROUPPMAP.(E).mLoc(:)');
            stdField = "LS";
            devField = "LD";
            stdLabel = "Local Standard";
            devLabel = "Local Deviant";
        else
            pCT = GROUPPMAP.(E).pGlo;
            m   = logical(GROUPPMAP.(E).mGlo(:)');
            stdField = "GS";
            devField = "GD";
            stdLabel = "Global Standard";
            devLabel = "Global Deviant";
        end

        t = double(GROUPPMAP.(E).times(:)');
        [rs, re] = supp2_find_runs(m);

        if isempty(rs)
            for ri = 1:numel(roiChans)
                roi = roiChans(ri);

                Troi = [Troi; table( ...
                    expLabels(ei), E, effectName, roi, false, ...
                    NaN, NaN, NaN, ...
                    NaN, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, ...
                    NaN, NaN, NaN, NaN, NaN, ...
                    stdLabel, devLabel, ...
                    'VariableNames', { ...
                    'ExperimentLabel','ExperimentKey','Effect','ROI','HasScalpWindow', ...
                    'WindowStartMS','WindowEndMS','WindowIndex', ...
                    'Nsub','StdMean_uV','StdSE_uV','DevMean_uV','DevSE_uV','DiffMean_uV','DiffSE_uV', ...
                    'tValue','df','pValue_Raw','Cohens_dz','MinPFWE_inROIWindow', ...
                    'StdLabel','DevLabel'})];
            end
            continue;
        end

        for w = 1:numel(rs)
            t1 = t(rs(w));
            t2 = t(re(w));
            idxWin = rs(w):re(w);

            for ri = 1:numel(roiChans)
                roi = roiChans(ri);
                roiField = matlab.lang.makeValidName(char(roi));

                if ~isfield(roiData, roiField) || isempty(roiData.(roiField).C)
                    continue;
                end

                Ccell = roiData.(roiField).C;
                nS = numel(Ccell);

                stdVals = nan(nS,1);
                devVals = nan(nS,1);

                for k = 1:nS
                    C = Ccell{k};
                    stdVals(k) = supp2_mean_in_window(C.relTimes, C.(stdField), t1, t2);
                    devVals(k) = supp2_mean_in_window(C.relTimes, C.(devField), t1, t2);
                end

                good = isfinite(stdVals) & isfinite(devVals);
                stdVals = stdVals(good);
                devVals = devVals(good);
                n = numel(stdVals);

                if n >= 2
                    [~, pRaw, ~, ST] = ttest(devVals, stdVals);

                    diffVals = devVals - stdVals;
                    diffMu = mean(diffVals, 'omitnan');
                    diffSE = std(diffVals, 0, 'omitnan') / sqrt(n);

                    if std(diffVals, 0, 'omitnan') > 0
                        dz = diffMu / std(diffVals, 0, 'omitnan');
                    else
                        dz = NaN;
                    end

                    tVal = ST.tstat;
                    df   = ST.df;
                else
                    pRaw = NaN;
                    tVal = NaN;
                    df   = NaN;
                    dz   = NaN;
                    diffMu = NaN;
                    diffSE = NaN;
                end

                chIdx = find(strcmpi(chanLabels, char(roi)), 1);
                if isempty(chIdx)
                    pMinROI = NaN;
                else
                    pMinROI = min(pCT(chIdx, idxWin), [], 'omitnan');
                end

                Troi = [Troi; table( ...
                    expLabels(ei), E, effectName, roi, true, ...
                    t1, t2, w, ...
                    n, ...
                    mean(stdVals, 'omitnan'), std(stdVals,0,'omitnan')/sqrt(max(n,1)), ...
                    mean(devVals, 'omitnan'), std(devVals,0,'omitnan')/sqrt(max(n,1)), ...
                    diffMu, diffSE, ...
                    tVal, df, pRaw, dz, pMinROI, ...
                    stdLabel, devLabel, ...
                    'VariableNames', { ...
                    'ExperimentLabel','ExperimentKey','Effect','ROI','HasScalpWindow', ...
                    'WindowStartMS','WindowEndMS','WindowIndex', ...
                    'Nsub','StdMean_uV','StdSE_uV','DevMean_uV','DevSE_uV','DiffMean_uV','DiffSE_uV', ...
                    'tValue','df','pValue_Raw','Cohens_dz','MinPFWE_inROIWindow', ...
                    'StdLabel','DevLabel'})];
            end
        end
    end
end

%% ------------------------------------------------------------
% S4: INDIVIDUAL-LEVEL SSI SUMMARY
%% ------------------------------------------------------------
Tindiv = table();

for ei = 1:numel(expList)
    E = expList(ei);
    if ~isfield(SSI, E)
        continue;
    end

    t = double(SSI.(E).times(:)');
    dt = median(diff(t));

    for fi = 1:numel(effList)
        effectName = effList(fi);

        if effectName == "Local"
            M = double(SSI.(E).local);
        else
            M = double(SSI.(E).global);
        end

        nSub = size(M,1);
        if nSub == 0
            continue;
        end

        hasEff = any(M > 0, 2);
        nWith  = sum(hasEff);
        pctWith = 100 * nWith / nSub;

        meanSSI = mean(M, 1, 'omitnan');
        [peakSSI, ixPeak] = max(meanSSI);
        tPeak = t(ixPeak);

        durMS = nan(nSub,1);
        for s = 1:nSub
            durMS(s) = sum(M(s,:) > 0) * dt;
        end

        if any(hasEff)
            medEffOnly = median(durMS(hasEff), 'omitnan');
        else
            medEffOnly = NaN;
        end

        Tindiv = [Tindiv; table( ...
            expLabels(ei), E, effectName, ...
            nSub, nWith, pctWith, ...
            peakSSI, tPeak, ...
            median(durMS, 'omitnan'), medEffOnly, ...
            'VariableNames', { ...
            'ExperimentLabel','ExperimentKey','Effect', ...
            'Nsub','Nsubjects_WithAnyEffect','PercentSubjects_WithAnyEffect', ...
            'PeakMeanSSI','PeakMeanSSI_TimeMS', ...
            'MedianDuration_AllSubjects_MS','MedianDuration_EffectOnly_MS'})];
    end
end

%% ------------------------------------------------------------
% S3 / S3b: SECONDARY Pz COMPARISON TABLES (from rainbox output)
%% ------------------------------------------------------------
Tsecondary = table();
TsecondarySubjects = table();

if ~isempty(OUT_EEG_RAINBOX_PZ) && isfield(OUT_EEG_RAINBOX_PZ, 'subjectTable') && isfield(OUT_EEG_RAINBOX_PZ, 'statsTable')
    Tsecondary = supp2_build_secondary_pz_table(OUT_EEG_RAINBOX_PZ);
    TsecondarySubjects = OUT_EEG_RAINBOX_PZ.subjectTable;
end

%% ------------------------------------------------------------
% S5: EEG INTERACTION SUMMARY
%% ------------------------------------------------------------
TinteractionEEG = supp2_build_interaction_table(INT, expLabels);


%% ------------------------------------------------------------
% SAVE CSVs
%% ------------------------------------------------------------
fileS1a = fullfile(tablesDir, 'EEG_S1a_Preprocessing_Long.csv');
fileS1b = fullfile(tablesDir, 'EEG_S1b_Preprocessing_Summary.csv');
fileS2  = fullfile(tablesDir, 'EEG_S2_Primary_GroupScalpEffects.csv');
fileS2b = fullfile(tablesDir, 'EEG_S2b_Primary_ROI_FollowUp.csv');
fileS3  = fullfile(tablesDir, 'EEG_S3_SecondaryPzComparisons.csv');
fileS3b = fullfile(tablesDir, 'EEG_S3b_SecondaryPz_SubjectValues.csv');
fileS4  = fullfile(tablesDir, 'EEG_S4_IndividualSSI.csv');
fileS5  = fullfile(tablesDir, 'EEG_S5_InteractionSummary.csv');

writetable(TretLong,    fileS1a);
writetable(TretSummary, fileS1b);
writetable(Tgroup,      fileS2);
writetable(Troi,        fileS2b);
if ~isempty(Tsecondary)
    writetable(Tsecondary, fileS3);
end
if ~isempty(TsecondarySubjects)
    writetable(TsecondarySubjects, fileS3b);
end
writetable(Tindiv,      fileS4);
writetable(TinteractionEEG, fileS5);


%% ------------------------------------------------------------
% SAVE ONE EXCEL WORKBOOK
%% ------------------------------------------------------------
xlsxFile = fullfile(tablesDir, 'EEG_Supplementary_Tables.xlsx');

writetable(TretLong,    xlsxFile, 'Sheet', 'S1a_Preproc_Long',   'FileType', 'spreadsheet');
writetable(TretSummary, xlsxFile, 'Sheet', 'S1b_Preproc_Sum',    'FileType', 'spreadsheet');
writetable(Tgroup,      xlsxFile, 'Sheet', 'S2_GroupScalp',      'FileType', 'spreadsheet');
writetable(Troi,        xlsxFile, 'Sheet', 'S2b_ROI',            'FileType', 'spreadsheet');
if ~isempty(Tsecondary)
    writetable(Tsecondary, xlsxFile, 'Sheet', 'S3_SecondaryStats', 'FileType', 'spreadsheet');
end
if ~isempty(TsecondarySubjects)
    writetable(TsecondarySubjects, xlsxFile, 'Sheet', 'S3b_SubjectVals', 'FileType', 'spreadsheet');
end
writetable(Tindiv,      xlsxFile, 'Sheet', 'S4_IndividualSSI',   'FileType', 'spreadsheet');
writetable(TinteractionEEG, xlsxFile, 'Sheet', 'S5_Interaction', 'FileType', 'spreadsheet');

%% ------------------------------------------------------------
% SAVE ONE WORD DOCUMENT
%% ------------------------------------------------------------
docxFile = fullfile(tablesDir, 'EEG_Supplementary_Tables.docx');

tablesForWord = struct();
tablesForWord(1).title = 'Supplementary Table S1a. EEG preprocessing / retention (long format)';
tablesForWord(1).T = TretLong;

tablesForWord(2).title = 'Supplementary Table S1b. EEG preprocessing / retention (summary)';
tablesForWord(2).T = TretSummary;

tablesForWord(3).title = 'Supplementary Table S2. Primary EEG group scalp-level effects';
tablesForWord(3).T = Tgroup;

tablesForWord(4).title = 'Supplementary Table S2b. Primary EEG ROI follow-up (Pz, Fz)';
tablesForWord(4).T = Troi;

if ~isempty(Tsecondary)
    tablesForWord(5).title = 'Supplementary Table S3. Secondary Pz between-task comparisons';
    tablesForWord(5).T = Tsecondary;
end

if ~isempty(TsecondarySubjects)
    tablesForWord(6).title = 'Supplementary Table S3b. Secondary Pz subject values';
    tablesForWord(6).T = TsecondarySubjects;
end

tablesForWord(7).title = 'Supplementary Table S4. Individual-level EEG SSI summary';
tablesForWord(7).T = Tindiv;

tablesForWord(8).title = 'Supplementary Table S5. EEG interaction summary';
tablesForWord(8).T = TinteractionEEG;

supp2_write_tables_to_word(tablesForWord, docxFile);

%% ------------------------------------------------------------
% OUTPUT
%% ------------------------------------------------------------
OUT = struct();
OUT.TretLong           = TretLong;
OUT.TretSummary        = TretSummary;
OUT.Tgroup             = Tgroup;
OUT.Troi               = Troi;
OUT.Tsecondary         = Tsecondary;
OUT.TsecondarySubjects = TsecondarySubjects;
OUT.Tindiv             = Tindiv;
OUT.Tinteraction       = TinteractionEEG;

OUT.files = struct();
OUT.files.tablesDir = tablesDir;
OUT.files.xlsx = xlsxFile;
OUT.files.docx = docxFile;
OUT.files.S1a_csv = fileS1a;
OUT.files.S1b_csv = fileS1b;
OUT.files.S2_csv  = fileS2;
OUT.files.S2b_csv = fileS2b;
OUT.files.S3_csv  = fileS3;
OUT.files.S3b_csv = fileS3b;
OUT.files.S4_csv  = fileS4;
OUT.files.S5_csv  = fileS5;

fprintf('\n[SUPPLEMENTARY EEG TABLES] Saved in:\n  %s\n', tablesDir);
fprintf('  Workbook: %s\n', xlsxFile);
fprintf('  Word:     %s\n', docxFile);

end


%% ============================================================
% HELPERS
%% ============================================================

function QA = supp2_load_info_table(outroot, grandFold)

QA = struct();
QA.exists = false;
QA.T = table();

cand = { ...
    fullfile(outroot,  'Info_EEG_LG.mat'), ...
    fullfile(grandFold,'Info_EEG_LG.mat')};

foundFile = '';
for i = 1:numel(cand)
    if exist(cand{i}, 'file')
        foundFile = cand{i};
        break;
    end
end

if isempty(foundFile)
    return;
end

S = load(foundFile);
fns = fieldnames(S);

if isfield(S, 'T_LG') && istable(S.T_LG)
    QA.T = S.T_LG;
    QA.exists = true;
    QA.file = foundFile;
    return;
end

for i = 1:numel(fns)
    if istable(S.(fns{i}))
        QA.T = S.(fns{i});
        QA.exists = true;
        QA.file = foundFile;
        return;
    end
end
end


function row = supp2_lookup_qa_row(QA, subjectID, expKey)

row = struct();
row.RejectedComponents = NaN;
row.RejectedTrials     = NaN;
row.GoodTrials         = NaN;
row.RejectedChannels   = NaN;
row.InfoSource         = "not_found";

if ~isfield(QA,'exists') || ~QA.exists || isempty(QA.T)
    return;
end

T = QA.T;
vn = string(T.Properties.VariableNames);

% ------------------------------------------------------------
% Step 1: try to match the subject row robustly
% ------------------------------------------------------------
rowIdx = [];

% subject number from 'S01' -> 1, 'S12' -> 12, etc.
tok = regexp(char(subjectID), '^S(\d+)$', 'tokens', 'once');
subNum = NaN;
if ~isempty(tok)
    subNum = str2double(tok{1});
end

% (a) first try obvious subject columns
candCols = ["Subject","subject","Sub","sub","Subject_ID","participant","Participant"];
for c = 1:numel(candCols)
    if ismember(candCols(c), vn)
        x = T.(candCols(c));

        xStr = string(x);

        m = strcmpi(strtrim(xStr), string(subjectID));
        if ~any(m)
            m = strcmpi(erase(strtrim(xStr)," "), erase(strtrim(string(subjectID))," "));
        end
        if ~any(m) && isfinite(subNum)
            xNum = suppress_num_convert(x);
            if ~all(isnan(xNum))
                m = (xNum == subNum);
            end
        end

        if any(m)
            rowIdx = find(m,1);
            break;
        end
    end
end

% (b) if still not found, scan ALL columns for the subject number or label
if isempty(rowIdx)
    for c = 1:numel(vn)
        x = T.(vn(c));

        % try string match
        try
            xStr = string(x);
            m = strcmpi(strtrim(xStr), string(subjectID));
            if ~any(m)
                m = strcmpi(erase(strtrim(xStr)," "), erase(strtrim(string(subjectID))," "));
            end
            if any(m)
                rowIdx = find(m,1);
                break;
            end
        catch
        end

        % try numeric subject number
        if isfinite(subNum)
            xNum = suppress_num_convert(x);
            if ~all(isnan(xNum))
                m = (xNum == subNum);
                if any(m)
                    rowIdx = find(m,1);
                    break;
                end
            end
        end
    end
end

if isempty(rowIdx)
    return;
end

r = T(rowIdx,:);

% ------------------------------------------------------------
% Step 2: direct field mapping
% ------------------------------------------------------------
badChVar = "";

if ismember("Num_Bad_Channels", vn)
    badChVar = "Num_Bad_Channels";
elseif ismember("Num_Bad_channels", vn)
    badChVar = "Num_Bad_channels";
elseif ismember("Num_bad_channels", vn)
    badChVar = "Num_bad_channels";
else
    idxBadCh = find(contains(lower(vn),"bad") & contains(lower(vn),"channel"), 1, 'first');
    if ~isempty(idxBadCh)
        badChVar = vn(idxBadCh);
    end
end

if badChVar ~= ""
    row.RejectedChannels = supp2_to_num(r.(badChVar));
end

if expKey == "E1" && ismember("Num_bad_trials_E1", vn)
    row.RejectedTrials = double(r.("Num_bad_trials_E1"));
elseif expKey == "E2" && ismember("Num_bad_trials_E2", vn)
    row.RejectedTrials = double(r.("Num_bad_trials_E2"));
end

if expKey == "E1" && ismember("Num_good_trials_E1", vn)
    row.GoodTrials = double(r.("Num_good_trials_E1"));
elseif expKey == "E2" && ismember("Num_good_trials_E2", vn)
    row.GoodTrials = double(r.("Num_good_trials_E2"));
end

if ismember("Bad_components_concatenated", vn)
    bc = r.("Bad_components_concatenated");

    if iscell(bc)
        bc = bc{1};
    end

    if isnumeric(bc) || islogical(bc)
        row.RejectedComponents = numel(bc);
    elseif isstring(bc) || ischar(bc)
        tok = regexp(char(bc), 'x(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            row.RejectedComponents = str2double(tok{1});
        end
    end
end

row.InfoSource = "Info_EEG_LG.mat";
end

function xNum = suppress_num_convert(x)
try
    if isnumeric(x)
        xNum = double(x);
        xNum = xNum(:);
    elseif iscell(x)
        xNum = nan(numel(x),1);
        for i = 1:numel(x)
            if isnumeric(x{i}) && isscalar(x{i})
                xNum(i) = double(x{i});
            elseif ischar(x{i}) || isstring(x{i})
                tmp = str2double(string(x{i}));
                if isfinite(tmp), xNum(i) = tmp; end
            end
        end
    else
        tmp = str2double(string(x));
        xNum = tmp(:);
    end
catch
    xNum = NaN(size(x,1),1);
end
end


function idx = supp2_find_var(vnl, candidates)
idx = [];
for i = 1:numel(candidates)
    j = find(contains(vnl, candidates{i}), 1, 'first');
    if ~isempty(j)
        idx = j;
        return;
    end
end
end


function idx = supp2_find_var_combo(vnl, mustContainOne, mustContainTwo)
idx = [];
for i = 1:numel(vnl)
    c1 = any(contains(vnl(i), mustContainOne));
    c2 = any(contains(vnl(i), mustContainTwo));
    if c1 && c2
        idx = i;
        return;
    end
end
end


function x = supp2_to_num(v)
try
    if isnumeric(v)
        x = double(v);
    elseif islogical(v)
        x = double(v);
    else
        x = str2double(string(v));
    end
catch
    x = NaN;
end
end


function subs = supp2_list_subjects_with_file(outroot, expKey, excludeSubs)
d = dir(fullfile(outroot, 'S*'));
d = d([d.isdir]);

names = string({d.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+$', 'once'));
subs = names(keep);

if ~isempty(excludeSubs)
    subs = subs(~ismember(lower(subs), lower(string(excludeSubs))));
end

hasFile = false(size(subs));
for i = 1:numel(subs)
    f = fullfile(outroot, char(subs(i)), sprintf('auditory_all_LG_%s_%s.set', char(subs(i)), char(expKey)));
    hasFile(i) = exist(f,'file') == 2;
end

subs = subs(hasFile);

nums = nan(size(subs));
for i = 1:numel(subs)
    tok = regexp(char(subs(i)), '^S(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        nums(i) = str2double(tok{1});
    end
end
[~,ix] = sort(nums);
subs = subs(ix);
end


function counts = supp2_count_trials_from_labels(EEG, targets)

counts = struct('FLS',NaN,'RLS',NaN,'FLD',NaN,'RLD',NaN, ...
    'LocalStandard',NaN,'LocalDeviant',NaN,'GlobalStandard',NaN,'GlobalDeviant',NaN);

if ~(isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type'))
    return;
end

M = false(4, EEG.trials);
for ep = 1:EEG.trials
    labs = supp2_to_strlist(EEG.epoch(ep).eventlg_type);
    for k = 1:4
        M(k,ep) = any(strcmpi(strtrim(labs), targets(k)));
    end
end

counts.FLS = sum(M(1,:));
counts.RLS = sum(M(2,:));
counts.FLD = sum(M(3,:));
counts.RLD = sum(M(4,:));

counts.LocalStandard  = counts.FLS + counts.RLS;
counts.LocalDeviant   = counts.FLD + counts.RLD;
counts.GlobalStandard = counts.FLS + counts.FLD;
counts.GlobalDeviant  = counts.RLS + counts.RLD;
end


function chanLabels = supp2_get_template_chanlabels(outroot, expKey, excludeSubs)
chanLabels = {};
subs = supp2_list_subjects_with_file(outroot, expKey, excludeSubs);

for i = 1:numel(subs)
    f = fullfile(outroot, char(subs(i)), sprintf('auditory_all_LG_%s_%s.set', char(subs(i)), char(expKey)));
    if exist(f,'file')
        EEG = pop_loadset(f);
        EEG = eeg_checkset(EEG);
        chanLabels = {EEG.chanlocs.labels};
        return;
    end
end
end


function C = supp2_load_subject_curves_eeg(setFile, chanName, targets)
C = [];

EEG = pop_loadset(setFile);
EEG = eeg_checkset(EEG);

if ~(isfield(EEG,'epoch') && isfield(EEG.epoch,'eventlg_type'))
    return;
end

ch = find(strcmpi({EEG.chanlocs.labels}, char(chanName)), 1);
if isempty(ch)
    return;
end

M = false(4, EEG.trials);
for ep = 1:EEG.trials
    labs = supp2_to_strlist(EEG.epoch(ep).eventlg_type);
    for k = 1:4
        M(k,ep) = any(strcmpi(strtrim(labs), targets(k)));
    end
end

idxFLS = find(M(1,:));
idxRLS = find(M(2,:));
idxFLD = find(M(3,:));
idxRLD = find(M(4,:));

if isempty(idxFLS) || isempty(idxRLS) || isempty(idxFLD) || isempty(idxRLD)
    return;
end

FLS = squeeze(mean(double(EEG.data(ch,:,idxFLS)),3)); FLS = FLS(:)';
RLS = squeeze(mean(double(EEG.data(ch,:,idxRLS)),3)); RLS = RLS(:)';
FLD = squeeze(mean(double(EEG.data(ch,:,idxFLD)),3)); FLD = FLD(:)';
RLD = squeeze(mean(double(EEG.data(ch,:,idxRLD)),3)); RLD = RLD(:)';

C = struct();
C.relTimes = double(EEG.times(:))';
C.FLS = FLS; C.RLS = RLS; C.FLD = FLD; C.RLD = RLD;
C.LS  = mean([FLS; RLS], 1, 'omitnan');
C.LD  = mean([FLD; RLD], 1, 'omitnan');
C.GS  = mean([FLS; FLD], 1, 'omitnan');
C.GD  = mean([RLS; RLD], 1, 'omitnan');
end


function S = supp2_to_strlist(v)
if isstring(v)
    S = v;
elseif ischar(v)
    S = string(v);
elseif iscell(v)
    S = strings(0,1);
    for i = 1:numel(v)
        if isstring(v{i}) || ischar(v{i})
            S(end+1,1) = string(v{i}); %#ok<AGROW>
        end
    end
else
    S = strings(0,1);
end
S = strtrim(regexprep(S,'\s+',' '));
end


function val = supp2_mean_in_window(relTimes, curve, startMS, endMS)
idx = relTimes >= startMS & relTimes <= endMS;
if ~any(idx) || all(isnan(curve(idx)))
    val = NaN;
else
    val = mean(curve(idx), 'omitnan');
end
end


function [rs, re] = supp2_find_runs(mask)
mask = logical(mask(:)');
d = diff([false mask false]);
rs = find(d == 1);
re = find(d == -1) - 1;
end


function Tsecondary = supp2_build_secondary_pz_table(OUTPZ)

Tsub   = OUTPZ.subjectTable;
Tstats = OUTPZ.statsTable;
cfgRB  = OUTPZ.cfg;

rows = {};
for effectName = ["Local","Global"]
    for measureName = ["Peak latency","Mean amplitude"]
        for condName = ["Standard","Deviant","Deviant - Standard contrast"]
            rows(end+1,:) = {effectName, measureName, condName}; %#ok<AGROW>
        end
    end
end

Tsecondary = table();

for i = 1:size(rows,1)
    effectName  = string(rows{i,1});
    measureName = string(rows{i,2});
    condSimple  = string(rows{i,3});

    [xE1, xE2, condFull, unitLabel] = supp2_get_secondary_vectors(Tsub, effectName, measureName, condSimple);
    good = isfinite(xE1) & isfinite(xE2);
    xE1 = xE1(good);
    xE2 = xE2(good);

    n = numel(xE1);

    if n >= 2
        diffVals = xE1 - xE2;
        diffMean = mean(diffVals, 'omitnan');
        diffSE   = std(diffVals,0,'omitnan') / sqrt(n);

        if std(diffVals,0,'omitnan') > 0
            dz = diffMean / std(diffVals,0,'omitnan');
        else
            dz = NaN;
        end
    else
        diffMean = NaN;
        diffSE   = NaN;
        dz       = NaN;
    end

    idx = (Tstats.Effect == effectName) & ...
        (Tstats.Measure == measureName) & ...
        (Tstats.Condition == condFull);

    if any(idx)
        rawP  = Tstats.RawPValue(find(idx,1,'first'));
        holmP = Tstats.CorrectedPValue(find(idx,1,'first'));
        sig   = Tstats.Significance(find(idx,1,'first'));
    else
        rawP  = NaN;
        holmP = NaN;
        sig   = "";
    end

    Tsecondary = [Tsecondary; table( ...
        effectName, measureName, condFull, ...
        string(cfgRB.statsMethod), n, ...
        mean(xE1,'omitnan'), std(xE1,0,'omitnan')/sqrt(max(n,1)), ...
        mean(xE2,'omitnan'), std(xE2,0,'omitnan')/sqrt(max(n,1)), ...
        diffMean, diffSE, dz, rawP, holmP, sig, unitLabel, ...
        'VariableNames', { ...
        'Effect','Measure','Condition','InferentialTest', ...
        'Npairs', ...
        'Counting_Mean','Counting_SE', ...
        'VisualInterference_Mean','VisualInterference_SE', ...
        'DifferenceMean_CountingMinusVI','DifferenceSE_CountingMinusVI','Cohens_dz', ...
        'RawPValue','HolmPValue','Significance','Unit'})];
end
end


function [xE1, xE2, condFull, unitLabel] = supp2_get_secondary_vectors(Tsub, effectName, measureName, condSimple)

switch effectName
    case "Local"
        prefix = "Local";
    otherwise
        prefix = "Global";
end

switch measureName
    case "Peak latency"
        suffixE1_std = prefix + "Standard_Latency_E1_ms";
        suffixE2_std = prefix + "Standard_Latency_E2_ms";
        suffixE1_dev = prefix + "Deviant_Latency_E1_ms";
        suffixE2_dev = prefix + "Deviant_Latency_E2_ms";
        unitLabel = "ms";
    otherwise
        suffixE1_std = prefix + "Standard_Mean_E1_uV";
        suffixE2_std = prefix + "Standard_Mean_E2_uV";
        suffixE1_dev = prefix + "Deviant_Mean_E1_uV";
        suffixE2_dev = prefix + "Deviant_Mean_E2_uV";
        unitLabel = "uV";
end

xStdE1 = Tsub.(suffixE1_std);
xStdE2 = Tsub.(suffixE2_std);
xDevE1 = Tsub.(suffixE1_dev);
xDevE2 = Tsub.(suffixE2_dev);

switch condSimple
    case "Standard"
        xE1 = xStdE1;
        xE2 = xStdE2;
        condFull = prefix + " Standard";
    case "Deviant"
        xE1 = xDevE1;
        xE2 = xDevE2;
        condFull = prefix + " Deviant";
    otherwise
        xE1 = xDevE1 - xStdE1;
        xE2 = xDevE2 - xStdE2;
        condFull = "Deviant - Standard contrast";
end
end


function supp2_write_tables_to_word(tablesForWord, docxFile)

try
    word = actxserver('Word.Application');
    word.Visible = false;
    doc = word.Documents.Add;
    sel = word.Selection;

    for i = 1:numel(tablesForWord)
        titleStr = tablesForWord(i).title;
        T = tablesForWord(i).T;

        sel.Font.Bold = 1;
        sel.Font.Size = 12;
        sel.TypeText(char(titleStr));
        sel.TypeParagraph;
        sel.TypeParagraph;

        nRows = height(T) + 1;
        nCols = width(T);

        rng = sel.Range;
        tbl = doc.Tables.Add(rng, nRows, nCols);

        % header
        for c = 1:nCols
            tbl.Cell(1,c).Range.Text = char(string(T.Properties.VariableNames{c}));
        end

        % body
        for r = 1:height(T)
            for c = 1:nCols
                v = T{r,c};
                tbl.Cell(r+1,c).Range.Text = supp2_val_to_text(v);
            end
        end

        tbl.Borders.Enable = 1;
        sel.MoveDown;
        sel.TypeParagraph;
        sel.TypeParagraph;
    end

    doc.SaveAs2(docxFile, 16); % wdFormatDocumentDefault = 16 (.docx)
    doc.Close(false);
    word.Quit;

catch
    warning('Word export failed.');
end
end


function txt = supp2_val_to_text(v)

try
    if iscell(v)
        if isempty(v)
            txt = "";
        else
            txt = supp2_val_to_text(v{1});
        end
        return;
    end

    if isstring(v)
        txt = char(join(v, ", "));
        return;
    end

    if ischar(v)
        txt = v;
        return;
    end

    if islogical(v)
        txt = char(string(v));
        return;
    end

    if isnumeric(v)
        if isempty(v) || any(isnan(v(:)))
            txt = "NaN";
        elseif isscalar(v)
            txt = num2str(v, '%.6g');
        else
            txt = char(strjoin(arrayfun(@(x) string(num2str(x,'%.6g')), v(:)', 'UniformOutput', false), ", "));
        end
        return;
    end

    txt = char(string(v));

catch
    txt = "";
end
end

function [mu,se] = supp2_mean_se(x)
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

function TinteractionEEG = supp2_build_interaction_table(INT, expLabels)

TinteractionEEG = table();

expKeys = ["E1","E2"];
for i = 1:numel(expKeys)
    E = expKeys(i);

    if ~isstruct(INT) || ~isfield(INT, E) || ~isfield(INT.(E), 'grey_omnibus')
        row = table(expLabels(i), E, false, NaN, NaN, NaN, NaN, "", ...
            'VariableNames', {'ExperimentLabel','ExperimentKey','HasSignificantWindow','StartMS','EndMS','MinPFWE','MinPTimeMS','MinPChannel'});
        TinteractionEEG = [TinteractionEEG; row];
        continue;
    end

    t = INT.(E).times(:)';
    p = INT.(E).grey_omnibus.p;
    m = logical(INT.(E).grey_omnibus.mask(:)');

    d = diff([false m false]);
    rs = find(d==1);
    re = find(d==-1)-1;

    if isempty(rs)
        row = table(expLabels(i), E, false, NaN, NaN, NaN, NaN, "", ...
            'VariableNames', {'ExperimentLabel','ExperimentKey','HasSignificantWindow','StartMS','EndMS','MinPFWE','MinPTimeMS','MinPChannel'});
        TinteractionEEG = [TinteractionEEG; row];
    else
        idxWin = rs(1):re(end);
        pWin = p(:,idxWin);
        [pmin, linIdx] = min(pWin(:));
        [chIdx, tIdxLocal] = ind2sub(size(pWin), linIdx);
        tMin = t(idxWin(tIdxLocal));

        row = table(expLabels(i), E, true, t(rs(1)), t(re(end)), pmin, tMin, "Ch"+string(chIdx), ...
            'VariableNames', {'ExperimentLabel','ExperimentKey','HasSignificantWindow','StartMS','EndMS','MinPFWE','MinPTimeMS','MinPChannel'});
        TinteractionEEG = [TinteractionEEG; row];
    end
end
end


%%