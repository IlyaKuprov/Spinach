% creates the Map with experimental parameters
%
% input parameters:
% freqMeas: measurement frequency (MW) in GHz
% field: main field for ENDOR exp. in G
% deltaField: EPR step size in G
% res_EN: resolution in ENDOR exp. in MHz
% range_EN: ENDOR exp. range in MHz
% t: pulse times in ns
%
% output parameters:
% Expt: the Map with the experimental parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function Expt=kehl_exp_create(freqMeas,field,deltaField,res_EN,range_EN,t)

    % Check consistency
    grumble(freqMeas,field,deltaField,res_EN,range_EN,t);
    Expt=containers.Map;

    % GHz to Hz
    Expt("FreqMeas")=freqMeas*10^9;

    % G to T
    Expt("Field")=field*10^-4;

    % G to T
    Expt("deltaField")=deltaField*10^-4;

    % MHz to Hz
    Expt("res_EN")=res_EN*10^6;
    Expt("range_EN")=range_EN*10^6;

    % ns to s
    Expt("t")=t*10^-9;

end

function grumble(freqMeas,field,deltaField,res_EN,range_EN,t)
if ~isnumeric(freqMeas)
    error('freqMeas must be numeric.');
end
if ~isnumeric(field)
    error('field must be numeric.');
end
if ~isnumeric(deltaField)
    error('deltaField must be numeric.');
end
if ~isnumeric(res_EN)
    error('res_EN must be numeric.');
end
if ~isnumeric(range_EN)
    error('range_EN must be numeric.');
end
if ~isnumeric(t)
    error('t must be numeric.');
end
end

