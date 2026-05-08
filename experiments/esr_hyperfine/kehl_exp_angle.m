% defines turning angle of RF pulse
%
% input parameters:
% expt_in: the Map with the experimental parameters
% value: turning angle of pulse in degree
%
% output parameters:
% Expt_out: updated Map with the experimental parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)
%


function Expt_out=kehl_exp_angle(expt_in,value)

    % Check consistency
    grumble(expt_in,value);
    Expt_out=expt_in;
    Expt_out("ang")=value;
end

function grumble(expt_in,value)
if ~isa(expt_in,'containers.Map')
    error('expt_in must be a containers.Map object.');
end
if ~isnumeric(value)
    error('value must be numeric.');
end
end

