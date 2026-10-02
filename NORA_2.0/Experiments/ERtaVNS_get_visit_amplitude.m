function amp = ERtaVNS_get_visit_amplitude(participantNumber,visitNumber,sessionCondition,sessionSite,askAlways)
% Return the visit-level stimulation command amplitude in the SAME convention
% as the lab's working Experiment_NORA script.
%
% CONTROL calls with askAlways=true and therefore asks the operator.
% E1/E2 call with askAlways=false:
%   1) try to reuse amplitude from the latest valid control MAT for the visit;
%   2) if unavailable, try an earlier E1/E2 MAT from that visit;
%   3) otherwise ask manually.
%
% This helper never creates a separate configuration MAT file.

    if nargin < 5
        askAlways = false;
    end

    participantID = sprintf('S%02d',participantNumber);
    rootDir = fileparts(mfilename('fullpath'));
    visitDir = fullfile(rootDir,'Data',participantID,sprintf('Visit%d',visitNumber));

    amp = NaN;

    if ~askAlways && exist(visitDir,'dir') == 7
        patterns = { ...
            sprintf('%sV%dC*.mat',participantID,visitNumber), ...
            sprintf('%sV%dE1*.mat',participantID,visitNumber), ...
            sprintf('%sV%dE2*.mat',participantID,visitNumber)};

        for p = 1:numel(patterns)
            d = dir(fullfile(visitDir,patterns{p}));
            d = d(~contains({d.name},'_INTERRUPTED'));
            if isempty(d), continue; end
            [~,ord] = sort([d.datenum],'descend');
            d = d(ord);

            for k = 1:numel(d)
                try
                    S = load(fullfile(visitDir,d(k).name),'experimentInfo');
                    if ~isfield(S,'experimentInfo'), continue; end
                    ei = S.experimentInfo;

                    candidate = NaN;
                    if isfield(ei,'stim') && isstruct(ei.stim) && isfield(ei.stim,'amp')
                        candidate = double(ei.stim.amp);
                    elseif isfield(ei,'hardwareStatus') && isstruct(ei.hardwareStatus) && ...
                            isfield(ei.hardwareStatus,'vnsAmplitude')
                        candidate = double(ei.hardwareStatus.vnsAmplitude);
                    end

                    if isfinite(candidate) && candidate > 0 && candidate <= 5
                        amp = candidate;
                        fprintf('Reusing visit amplitude %.6g from %s\n',amp,d(k).name);
                        return;
                    end
                catch
                end
            end
        end
    end

    label = sprintf(['%s Visit %d | %s | %s\n' ...
        'Enter stimulation command amplitude in the SAME 0-5 convention used by Experiment_NORA: '], ...
        participantID,visitNumber,upper(char(sessionCondition)),char(sessionSite));

    raw = input(label,'s');
    if contains(raw,',')
        error('Use a dot (.) as the decimal separator, not a comma (,).');
    end
    amp = str2double(raw);

    if ~isfinite(amp) || amp <= 0 || amp > 5
        error('Stimulation command amplitude must be >0 and <=5 in the Experiment_NORA convention.');
    end
end
