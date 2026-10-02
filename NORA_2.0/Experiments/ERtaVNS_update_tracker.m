function ERtaVNS_update_tracker(participantNumber,visitNumber,component,status,filePath)
% Update the advisory participant tracker after a completed/interrupted run.
% Failure to update the tracker is NON-FATAL and never invalidates the saved
% experiment data.

    try
        if ~ismember(visitNumber,[1 2])
            warning('Tracker not updated: visit number must be 1 or 2.');
            return;
        end

        participantID = sprintf('S%02d',participantNumber);
        component = upper(strtrim(char(component)));
        status = upper(strtrim(char(status)));

        validComponents = {'CONTROL','E1','E2'};
        if ~ismember(component,validComponents)
            warning('Tracker not updated: unknown component %s.',component);
            return;
        end

        trackerFile = fullfile(fileparts(mfilename('fullpath')),'ERtaVNS_participant_tracker.csv');
        T = readtable(trackerFile,'TextType','string','VariableNamingRule','preserve');
        idx = find(upper(strtrim(string(T.Participant))) == participantID,1);

        if isempty(idx)
            warning('Tracker not updated: participant %s not found.',participantID);
            return;
        end

        if strcmp(component,'CONTROL')
            stem = sprintf('V%d_Control',visitNumber);
        else
            stem = sprintf('V%d_%s',visitNumber,component);
        end

        statusCol = [stem '_Status'];
        timeCol   = [stem '_Time'];
        fileCol   = [stem '_File'];

        T.(statusCol) = string(T.(statusCol));
        T.(timeCol)   = string(T.(timeCol));
        T.(fileCol)   = string(T.(fileCol));

        T.(statusCol)(idx) = string(status);
        T.(timeCol)(idx) = string(datestr(now,'yyyy-mm-dd HH:MM:SS'));
        T.(fileCol)(idx) = string(filePath);

        writetable(T,trackerFile);
        fprintf('Tracker updated: %s Visit %d %s = %s\n',participantID,visitNumber,component,status);

    catch ME
        txt = 'tracker update error';
        try
            txt = getReport(ME,'basic','hyperlinks','off');
            txt = regexprep(txt,'[\r\n\t]+',' ');
            txt = regexprep(txt,'\s+',' ');
        catch
        end
        warning('Experiment data were saved, but the tracker could not be updated: %s',strtrim(txt));
    end
end
