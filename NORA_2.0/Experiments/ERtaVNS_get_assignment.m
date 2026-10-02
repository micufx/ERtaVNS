function assignment = ERtaVNS_get_assignment(participantNumber,visitNumber)
% Read the pre-generated participant/session assignment from the tracker.
% The tracker is the source of truth for visit condition, task order and
% local-global physical block order.

    if ~ismember(visitNumber,[1 2])
        error('Visit number must be 1 or 2.');
    end

    participantID = sprintf('S%02d',participantNumber);
    trackerFile = fullfile(fileparts(mfilename('fullpath')),'ERtaVNS_participant_tracker.csv');

    if exist(trackerFile,'file') ~= 2
        error('Participant tracker not found: %s',trackerFile);
    end

    T = readtable(trackerFile,'TextType','string','VariableNamingRule','preserve');
    ids = upper(strtrim(string(T.Participant)));
    idx = find(ids == participantID,1);

    if isempty(idx)
        % Do not block acquisition merely because the requested positive ID is
        % beyond the pre-generated tracker. Use a deterministic 16-cell
        % counterbalancing cycle for new IDs. Existing tracker assignments are
        % never changed.
        scheduleGroup = mod(participantNumber-1,4)+1;
        blockGroup = mod((participantNumber-1)+floor((participantNumber-1)/4),4)+1;
        switch scheduleGroup
            case 1, v1Cond='ACTIVE'; v1Order={'E1','E2'}; v2Cond='SHAM';   v2Order={'E2','E1'};
            case 2, v1Cond='ACTIVE'; v1Order={'E2','E1'}; v2Cond='SHAM';   v2Order={'E1','E2'};
            case 3, v1Cond='SHAM';   v1Order={'E1','E2'}; v2Cond='ACTIVE'; v2Order={'E2','E1'};
            otherwise, v1Cond='SHAM';v1Order={'E2','E1'}; v2Cond='ACTIVE'; v2Order={'E1','E2'};
        end
        blockMap={{'a','d'},{'d','a'},{'b','c'},{'c','b'}};
        if visitNumber==1, cond=v1Cond; taskOrder=v1Order; else, cond=v2Cond; taskOrder=v2Order; end
        assignment=struct('participantID',participantID,'participantNumber',participantNumber, ...
            'visitNumber',visitNumber,'sessionCondition',cond,'taskOrder',{taskOrder}, ...
            'blockOrder',{blockMap{blockGroup}},'scheduleGroup',scheduleGroup, ...
            'blockGroup',blockGroup,'trackerRow',NaN);
        warning(['%s is beyond/not present in the pre-generated tracker. A deterministic ' ...
            'counterbalanced assignment was generated so acquisition can proceed.'],participantID);
        return;
    end

    if visitNumber == 1
        cond = upper(strtrim(string(T.Visit1Condition(idx))));
        orderText = upper(strtrim(string(T.Visit1TaskOrder(idx))));
    else
        cond = upper(strtrim(string(T.Visit2Condition(idx))));
        orderText = upper(strtrim(string(T.Visit2TaskOrder(idx))));
    end

    taskOrder = regexp(char(orderText),'>','split');
    taskOrder = taskOrder(~cellfun('isempty',taskOrder));
    if numel(taskOrder) ~= 2
        error('Invalid task order for %s Visit %d in tracker.',participantID,visitNumber);
    end

    blockText = lower(strtrim(string(T.BlockOrder(idx))));
    blockOrder = regexp(char(blockText),'>','split');
    blockOrder = blockOrder(~cellfun('isempty',blockOrder));
    if numel(blockOrder) ~= 2
        error('Invalid block order for %s in tracker.',participantID);
    end

    assignment = struct();
    assignment.participantID = participantID;
    assignment.participantNumber = participantNumber;
    assignment.visitNumber = visitNumber;
    assignment.sessionCondition = char(cond);
    assignment.taskOrder = taskOrder;
    assignment.blockOrder = blockOrder;
    assignment.scheduleGroup = double(T.ScheduleGroup(idx));
    assignment.blockGroup = double(T.BlockGroup(idx));
    assignment.trackerRow = idx;
end
