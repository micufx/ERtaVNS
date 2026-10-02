function T_LG = ensure_pupil_timing_columns(T_LG)

% Ensure_pupil_timing_columns
% Adds (if missing) the timing/QC columns used by the pupil preprocessing script.
%
% Columns created / ensured:
%   SOA_E1_Pupil, SOA_E2_Pupil                 (double)
%   ITI_E1_Pupil, ITI_E2_Pupil                 (double)
%   BaselineMapped_E1_Pupil, BaselineMapped_E2_Pupil   (cell) each row: 1x2 double
%   EpochSuggest_E1_Pupil,   EpochSuggest_E2_Pupil     (cell) each row: 1x2 double
%   BaselineAnchor_E1_Pupil, BaselineAnchor_E2_Pupil   (string)
%
% NOTE:
% Cell columns are used so the table displays numeric pairs like:
%   [-2179 -1679]
% without quotes (unlike strings which show "[-2179 -1679]").

%% numeric columns (unchanged)
needNum = {'SOA_E1_Pupil','SOA_E2_Pupil','ITI_E1_Pupil','ITI_E2_Pupil'};
for k = 1:numel(needNum)
    vn = needNum{k};
    if ~ismember(vn, T_LG.Properties.VariableNames)
        T_LG.(vn) = nan(height(T_LG),1);
    end
end

%% cell columns (1x2 doubles)
needCell = {'BaselineMapped_E1_Pupil','BaselineMapped_E2_Pupil', ...
            'EpochSuggest_E1_Pupil','EpochSuggest_E2_Pupil'};

for k = 1:numel(needCell)
    vn = needCell{k};

    if ~ismember(vn, T_LG.Properties.VariableNames)
        % create fresh cell column
        T_LG.(vn) = cell(height(T_LG),1);
        continue;
    end

    % If it exists but is NOT cell (e.g., string from older run), try to convert
    if ~iscell(T_LG.(vn))
        old = T_LG.(vn);

        % Convert string/char/cellstr like "[-2164 -1664]" into {[-2164 -1664]}
        if isstring(old) || iscellstr(old)
            oldS = string(old);
            newC = cell(height(T_LG),1);

            for r = 1:height(T_LG)
                s = strtrim(oldS(r));
                if s=="" || s=="<missing>"
                    newC{r} = [];
                else
                    % Remove brackets if present, then parse numbers
                    s2 = erase(erase(s,'['),']');
                    nums = sscanf(s2, '%f');
                    if numel(nums) >= 2
                        newC{r} = nums(1:2).';
                    else
                        newC{r} = [];
                    end
                end
            end

            T_LG.(vn) = newC;

        else
            % Unknown type -> safest is to replace with empty cell column
            tmp = cell(height(T_LG),1);
            T_LG.(vn) = tmp;
        end
    end
end

%% anchor columns (strings)
needStr = {'BaselineAnchor_E1_Pupil','BaselineAnchor_E2_Pupil'};
for k = 1:numel(needStr)
    vn = needStr{k};
    if ~ismember(vn, T_LG.Properties.VariableNames)
        T_LG.(vn) = strings(height(T_LG),1);
    end
end

end
