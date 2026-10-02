function [participantNumber,participantID,visitNumber] = ERtaVNS_resolve_participant(participantNumber,visitNumber,component,scriptDir)
% Task-aware participant/visit resolver.
% A folder by itself never counts as a recording. A component counts as
% recorded only when VisitN contains at least one FILE whose name begins:
%   SxxVnC, SxxVnE1, or SxxVnE2 (extension/type does not matter).
% Existing data never blocks acquisition: the operator can always keep the
% requested participant/visit and chooseAvailableRunStem() protects files.

component=upper(char(component));
if ~ismember(component,{'C','E1','E2'})
    error('Unknown component "%s". Expected C, E1, or E2.',component);
end
participantID=sprintf('S%02d',participantNumber);
dataRoot=fullfile(scriptDir,'Data');

s1=visitStatus(dataRoot,participantID,1);
s2=visitStatus(dataRoot,participantID,2);
cur = s1; other=s2;
if visitNumber==2, cur=s2; other=s1; end

% If THIS component has not been recorded in the requested visit, proceed
% silently even when the participant folder or other task files already exist.
if ~cur.(component)
    return;
end

componentName=prettyComponent(component);
statusText=sprintf(['%s already appears to contain %s data.\n\n' ...
    'Visit 1: Control=%s | Experiment 1=%s | Experiment 2=%s\n' ...
    'Visit 2: Control=%s | Experiment 1=%s | Experiment 2=%s\n\n'], ...
    sprintf('%s Visit %d',participantID,visitNumber),componentName, ...
    yn(s1.C),yn(s1.E1),yn(s1.E2),yn(s2.C),yn(s2.E1),yn(s2.E2));

% Suggestion hierarchy:
% 1) same component in the contrary visit if not yet recorded;
% 2) if both visits are fully complete, next participant;
% 3) otherwise keep current visit (other components remain to be run with
%    their corresponding scripts), while clearly reporting what is missing.
if ~other.(component)
    suggestedVisit=3-visitNumber;
    msg=[statusText sprintf(['The same component is not yet recorded in Visit %d.\n' ...
        'Suggested action: use Visit %d.\n\nYou may still keep Visit %d for an intentional rerun.'], ...
        suggestedVisit,suggestedVisit,visitNumber)];
    choice=questdlg(msg,'Existing component recording', ...
        sprintf('Use Visit %d',suggestedVisit),sprintf('Keep Visit %d',visitNumber),'Cancel',sprintf('Use Visit %d',suggestedVisit));
    if isempty(choice) || strcmp(choice,'Cancel')
        error('Experiment cancelled at participant/visit selection.');
    elseif strcmp(choice,sprintf('Use Visit %d',suggestedVisit))
        visitNumber=suggestedVisit;
        fprintf('%s already recorded for %s Visit %d. Switched to Visit %d as suggested.\n', ...
            componentName,participantID,3-suggestedVisit,suggestedVisit);
    else
        fprintf('Keeping %s Visit %d %s as requested. Existing files will not be overwritten.\n', ...
            participantID,visitNumber,componentName);
    end
    return;
end

if s1.complete && s2.complete
    nextNumber=participantNumber+1;
    while participantComplete(dataRoot,nextNumber)
        nextNumber=nextNumber+1;
    end
    suggestedID=sprintf('S%02d',nextNumber);
    msg=[statusText sprintf(['All three components are present in BOTH visits for %s.\n' ...
        'Suggested action: continue with next participant %s.\n\n' ...
        'You may still keep %s Visit %d for an intentional rerun.'], ...
        participantID,suggestedID,participantID,visitNumber)];
    choice=questdlg(msg,'Participant complete', ...
        ['Use ' suggestedID],sprintf('Keep %s V%d',participantID,visitNumber),'Cancel',['Use ' suggestedID]);
    if isempty(choice) || strcmp(choice,'Cancel')
        error('Experiment cancelled at participant selection.');
    elseif strcmp(choice,['Use ' suggestedID])
        participantNumber=nextNumber;
        participantID=suggestedID;
        visitNumber=1;
        fprintf('Both visits complete. Changed to suggested next participant: %s Visit 1.\n',participantID);
    else
        fprintf('Keeping %s Visit %d as requested for an intentional rerun. Existing files will not be overwritten.\n',participantID,visitNumber);
    end
    return;
end

missing1=missingText(s1); missing2=missingText(s2);
msg=[statusText sprintf(['This %s component exists in both visits, but the participant is not fully complete.\n' ...
    'Still missing: Visit 1: %s | Visit 2: %s\n\n' ...
    'Run the corresponding missing component script(s). You may keep this selection for an intentional rerun.'], ...
    componentName,missing1,missing2)];
choice=questdlg(msg,'Component already recorded',sprintf('Keep %s V%d',participantID,visitNumber),'Cancel',sprintf('Keep %s V%d',participantID,visitNumber));
if isempty(choice) || strcmp(choice,'Cancel')
    error('Experiment cancelled at participant/visit selection.');
end
fprintf('Keeping %s Visit %d %s as requested. Existing files will not be overwritten.\n',participantID,visitNumber,componentName);
end

function s=visitStatus(dataRoot,pid,v)
visitDir=fullfile(dataRoot,pid,sprintf('Visit%d',v));
s.C = hasPrefixFile(visitDir,sprintf('%sV%dC',pid,v));
s.E1= hasPrefixFile(visitDir,sprintf('%sV%dE1',pid,v));
s.E2= hasPrefixFile(visitDir,sprintf('%sV%dE2',pid,v));
s.complete=s.C && s.E1 && s.E2;
end

function tf=hasPrefixFile(folder,prefix)
tf=false;
if exist(folder,'dir')~=7, return; end
d=dir(fullfile(folder,[prefix '*']));
for k=1:numel(d)
    if ~d(k).isdir
        tf=true; return;
    end
end
end

function tf=participantComplete(dataRoot,n)
pid=sprintf('S%02d',n);
a=visitStatus(dataRoot,pid,1); b=visitStatus(dataRoot,pid,2);
tf=a.complete && b.complete;
end

function t=yn(tf)
if tf, t='RECORDED'; else, t='missing'; end
end

function t=prettyComponent(c)
switch c
    case 'C', t='Control';
    case 'E1', t='Experiment 1';
    case 'E2', t='Experiment 2';
end
end

function t=missingText(s)
a={};
if ~s.C, a{end+1}='Control'; end %#ok<AGROW>
if ~s.E1, a{end+1}='Experiment 1'; end %#ok<AGROW>
if ~s.E2, a{end+1}='Experiment 2'; end %#ok<AGROW>
if isempty(a), t='none'; else, t=strjoin(a,', '); end
end
