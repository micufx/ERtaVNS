function write_TLG_excel_exact(T_LG, excelFile)
% write_TLG_excel_exact
% Writes T_LG to Excel WITHOUT "column explosion" from cells containing vectors.
% Strategy:
%   - Create an export copy
%   - Convert any cell/array-like content into single scalar strings per row
%   - Write the resulting table

T_export = T_LG;

vn = T_export.Properties.VariableNames;

for c = 1:numel(vn)
    col = T_export.(vn{c});

    % ---------- CASE 1: cell column (most problematic) ----------
    if iscell(col)
        S = strings(height(T_export),1);

        for r = 1:height(T_export)
            v = col{r};

            if isempty(v)
                S(r) = "";
                continue;
            end

            % numeric vector/matrix inside a cell -> one bracket string
            if isnumeric(v)
                if isscalar(v)
                    S(r) = sprintf('%g', v);
                else
                    % flatten row-major
                    vv = v(:).';
                    S(r) = "[" + strjoin(string(vv), " ") + "]";
                end
                continue;
            end

            % char -> string
            if ischar(v)
                S(r) = string(v);
                continue;
            end

            % string scalar
            if isstring(v) && isscalar(v)
                S(r) = v;
                continue;
            end

            % string array -> join
            if isstring(v) && ~isscalar(v)
                S(r) = strjoin(v(:).', ", ");
                continue;
            end

            % cellstr -> join
            if iscellstr(v)
                S(r) = strjoin(string(v(:).'), ", ");
                continue;
            end

            % fallback: try converting to string safely
            try
                S(r) = string(v);
            catch
                S(r) = "";
            end
        end

        T_export.(vn{c}) = S;
        continue;
    end

    % ---------- CASE 2: string arrays (non-scalar) ----------
    if isstring(col) && ~isvector(col)
        % unlikely in tables, but just in case
        col = col(:);
    end

    % If a column is a string array where each row might be non-scalar, join rows
    if isstring(col) && height(T_export) == numel(col)
        % OK: string scalar per row already
        continue;
    end

    % ---------- CASE 3: categorical arrays (safe) ----------
    % categorical writes fine; no action needed.

    % ---------- CASE 4: numeric matrices in a variable (also problematic) ----------
    % If someone stored NxK numeric data directly as a table variable, Excel will
    % write K columns. If you want to prevent that, convert to bracket-string.
    if isnumeric(col) && ~isvector(col) && size(col,1) == height(T_export)
        % Each row is a vector -> string-ify rowwise
        S = strings(height(T_export),1);
        for r = 1:height(T_export)
            vv = col(r,:);
            S(r) = "[" + strjoin(string(vv), " ") + "]";
        end
        T_export.(vn{c}) = S;
    end
end

% Write the export table
writetable(T_export, excelFile);
end
