% sets RF start frequency
%
% input parameters:
% expt_in: the Map with the experimental parameters
% value: RF start frequency in MHz
%
% output parameters:
% Expt_out: updated Map with the experimental parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function Expt_out=kehl_exp_rfstart(expt_in,value)

    % Check consistency
    grumble(expt_in,value);
    Expt_out=expt_in;

    Expt_out("RF_start")=value*1e6;
end

function grumble(expt_in,value)
if ~isa(expt_in,'containers.Map')
    error('expt_in must be a containers.Map object.');
end
if ~isnumeric(value)
    error('value must be numeric.');
end
end

