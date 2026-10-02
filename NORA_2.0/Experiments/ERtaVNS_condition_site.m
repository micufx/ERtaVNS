function site = ERtaVNS_condition_site(sessionCondition)
% Physical electrode site for the visit-level condition.

    c = upper(strtrim(char(sessionCondition)));
    switch c
        case 'ACTIVE'
            site = 'left cymba conchae';
        case 'SHAM'
            site = 'left earlobe';
        otherwise
            error('Unknown session condition: %s',c);
    end
end
