function [nextTask,detail] = ERtaVNS_next_task(participantNumber,visitNumber)
% Advisory-only tracker helper. It NEVER blocks an acquisition script.
% Returns CONTROL first if not marked DONE, then the two tasks in the
% participant's assigned order.

    nextTask = '';
    detail = '';

    try
        assignment = ERtaVNS_get_assignment(participantNumber,visitNumber);
        participantID = assignment.participantID;
        trackerFile = fullfile(fileparts(mfilename('fullpath')),'ERtaVNS_participant_tracker.csv');
        T = readtable(trackerFile,'TextType','string','VariableNamingRule','preserve');

        idx = find(upper(strtrim(string(T.Participant))) == participantID,1);
        prefix = sprintf('V%d_',visitNumber);

        controlCol = [prefix 'Control_Status'];
        controlVals = T.(controlCol);
        controlStatus = upper(strtrim(string(controlVals(idx))));
        if controlStatus ~= "DONE"
            nextTask = 'CONTROL';
            detail = 'Control phase is not marked DONE.';
            return;
        end

        for k = 1:numel(assignment.taskOrder)
            task = upper(assignment.taskOrder{k});
            col = [prefix task '_Status'];
            vals = T.(col);
            status = upper(strtrim(string(vals(idx))));
            if status ~= "DONE"
                nextTask = task;
                detail = sprintf('%s is the next not-DONE component in the assigned order.',task);
                return;
            end
        end

        nextTask = '';
        detail = 'All assigned components are marked DONE.';
    catch
        nextTask = '';
        detail = 'Tracker status could not be read; acquisition remains independently runnable.';
    end
end
