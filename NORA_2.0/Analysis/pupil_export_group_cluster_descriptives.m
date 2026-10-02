function OUT = pupil_export_group_cluster_descriptives(outpath, grandOutFolder, STATS, excludeSubs, useField)

if nargin < 4 || isempty(excludeSubs)
    excludeSubs = strings(0,1);
end
if nargin < 5 || isempty(useField)
    useField = "PUPIL_Z";
end

excludeSubs = string(excludeSubs);
useField    = string(useField);

tablesDir = fullfile(grandOutFolder, 'tables');
if ~exist(tablesDir, 'dir')
    mkdir(tablesDir);
end

expList   = ["E1","E2"];
expLabels = ["Counting","Visual interference"];

Tall = table();

for e = 1:numel(expList)
    expKey = expList(e);
    expLab = expLabels(e);

    subNames = pgcd_list_subjects_with_file(outpath, expKey, excludeSubs);
    if isempty(subNames)
        continue;
    end

    for eff = ["Local","Global"]

        if eff == "Local"
            mask = logical(STATS.(expKey).maskLocal(:)');
            clustersField = "clustersLocal";
        else
            mask = logical(STATS.(expKey).maskGlobal(:)');
            clustersField = "clustersGlobal";
        end

        t = double(STATS.(expKey).relTimes(:)');

        if ~any(mask)
            newRow = table( ...
                expLab, eff, 0, NaN, NaN, NaN, ...
                NaN, NaN, NaN, NaN, ...
                NaN, NaN, NaN, NaN, NaN, NaN, ...
                NaN, NaN, NaN, NaN, NaN, ...
                'VariableNames', { ...
                'Experiment','Effect','Nsub','StartMS','EndMS','Best_pFWE', ...
                'Std_Mean','Std_SE','Dev_Mean','Dev_SE', ...
                'Contrast_Mean','Contrast_SE','CI_Low','CI_High','t_value','df', ...
                'p_uncorrected','dz','PeakDiff','PeakMS','tmax','tmaxMS','p_at_tmax'});
            Tall = [Tall; newRow];
            continue;
        end

        startMS = t(find(mask,1,'first'));
        endMS   = t(find(mask,1,'last'));

        pBest = NaN;
        if isfield(STATS.(expKey), clustersField)
            C = STATS.(expKey).(clustersField);
            if ~isempty(C)
                try
                    pBest = min([C.pFWE]);
                catch
                    pBest = NaN;
                end
            end
        end

        Xstd = [];
        Xdev = [];
        relTimesRef = [];

        for s = 1:numel(subNames)
            sID = subNames(s);
            f = fullfile(outpath, char(sID), sprintf('%s_%s_pupil_preproc.mat', sID, expKey));
            if ~exist(f,'file')
                continue;
            end

            S = pgcd_load_subject_curves(f, useField);
            if isempty(S)
                continue;
            end

            if isempty(relTimesRef)
                relTimesRef = S.relTimes(:)';
            end

            if eff == "Local"
                Xstd(end+1,:) = S.LS; %#ok<AGROW>
                Xdev(end+1,:) = S.LD; %#ok<AGROW>
            else
                Xstd(end+1,:) = S.GS; %#ok<AGROW>
                Xdev(end+1,:) = S.GD; %#ok<AGROW>
            end
        end

        if isempty(Xstd) || isempty(Xdev)
            continue;
        end

        idx = relTimesRef >= startMS & relTimesRef <= endMS;
        tw  = relTimesRef(idx);

        subjStd = mean(Xstd(:,idx), 2, 'omitnan');
        subjDev = mean(Xdev(:,idx), 2, 'omitnan');
        subjCon = subjDev - subjStd;

        n  = numel(subjCon);
        df = n - 1;

        muStd = mean(subjStd, 'omitnan');
        muDev = mean(subjDev, 'omitnan');
        muCon = mean(subjCon, 'omitnan');

        seStd = std(subjStd, 0, 'omitnan') / sqrt(n);
        seDev = std(subjDev, 0, 'omitnan') / sqrt(n);
        seCon = std(subjCon, 0, 'omitnan') / sqrt(n);

        [~, punc, ci, stats] = ttest(subjCon, 0);
        tval = stats.tstat;
        dz   = muCon / std(subjCon, 0, 'omitnan');

        grandDiff = mean(Xdev - Xstd, 1, 'omitnan');
        [~, ixPeak] = max(abs(grandDiff(idx)));
        peakDiff = grandDiff(find(idx,1,'first') + ixPeak - 1);
        peakMS   = tw(ixPeak);

        Dwin  = Xdev(:,idx) - Xstd(:,idx);
        tvals = mean(Dwin,1,'omitnan') ./ (std(Dwin,0,1,'omitnan') ./ sqrt(n));
        [~, ixTmax] = max(abs(tvals));
        tmax    = tvals(ixTmax);
        tmaxMS  = tw(ixTmax);
        pAtTmax = 2 * tcdf(-abs(tmax), df);

        newRow = table( ...
            expLab, eff, n, startMS, endMS, pBest, ...
            muStd, seStd, muDev, seDev, ...
            muCon, seCon, ci(1), ci(2), tval, df, ...
            punc, dz, peakDiff, peakMS, tmax, tmaxMS, pAtTmax, ...
            'VariableNames', { ...
            'Experiment','Effect','Nsub','StartMS','EndMS','Best_pFWE', ...
            'Std_Mean','Std_SE','Dev_Mean','Dev_SE', ...
            'Contrast_Mean','Contrast_SE','CI_Low','CI_High','t_value','df', ...
            'p_uncorrected','dz','PeakDiff','PeakMS','tmax','tmaxMS','p_at_tmax'});

        Tall = [Tall; newRow];
    end
end

csvFile  = fullfile(tablesDir, 'PUPIL_Table_S2b_Group_Cluster_Descriptives.csv');
xlsxFile = fullfile(tablesDir, 'PUPIL_Table_S2b_Group_Cluster_Descriptives.xlsx');

writetable(Tall, csvFile);
writetable(Tall, xlsxFile, 'Sheet', 'GroupClusterDescriptives');

OUT = struct();
OUT.table = Tall;
OUT.csv   = csvFile;
OUT.xlsx  = xlsxFile;

disp(Tall);
fprintf('\n[PUPIL GROUP DESCRIPTIVES] Saved:\n  %s\n  %s\n', csvFile, xlsxFile);

end

% -------------------------------------------------------------------------
function subs = pgcd_list_subjects_with_file(outpath, expKey, excludeSubs)
d = dir(fullfile(outpath, 'S*'));
d = d([d.isdir]);

names = string({d.name});
names = names(~ismember(names, [".","..","GRAND_AVG"]));
keep = ~cellfun('isempty', regexp(cellstr(names), '^S\d+$', 'once'));
subs = names(keep);

if ~isempty(excludeSubs)
    subs = subs(~ismember(subs, excludeSubs));
end

hasFile = false(size(subs));
for i = 1:numel(subs)
    f = fullfile(outpath, char(subs(i)), sprintf('%s_%s_pupil_preproc.mat', subs(i), expKey));
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

% -------------------------------------------------------------------------
function C = pgcd_load_subject_curves(matFile, useField)
C = [];

S = load(matFile, 'dataTableClean', 'relTimes');
if ~isfield(S,'dataTableClean') || isempty(S.dataTableClean) || height(S.dataTableClean)==0
    return;
end
if ~isfield(S,'relTimes') || isempty(S.relTimes)
    return;
end

T = S.dataTableClean;
relTimes = S.relTimes(:)';

C = struct();
C.relTimes = relTimes;

C.LS = pgcd_subject_mean_curve(T, useField, 'sound', "Local Standard");
C.LD = pgcd_subject_mean_curve(T, useField, 'sound', "Local Deviant");

typeVec = string(T.type);
effectVec = strings(height(T),1);
effectVec(startsWith(typeVec, "Frequent")) = "Global Standard";
effectVec(startsWith(typeVec, "Rare"))     = "Global Deviant";

T2 = T;
T2.effect = effectVec;

C.GS = pgcd_subject_mean_curve(T2, useField, 'effect', "Global Standard");
C.GD = pgcd_subject_mean_curve(T2, useField, 'effect', "Global Deviant");
end

% -------------------------------------------------------------------------
function m = pgcd_subject_mean_curve(T, fieldName, condVar, condValue)
fieldName = char(fieldName);
condVec = string(T.(condVar));
idx = (condVec == string(condValue));

if ~any(idx)
    if iscell(T.(fieldName)) && ~isempty(T.(fieldName))
        m = nan(1, numel(T.(fieldName){1}));
    else
        m = nan(1,1);
    end
    return;
end

trialCells = T.(fieldName)(idx);
X = vertcat(trialCells{:});
m = mean(X, 1, 'omitnan');
end