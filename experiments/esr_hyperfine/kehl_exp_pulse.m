% sets pulse profile input file
%
% input parameters:
% expt_in: the Map with the experimental parameters
% data: filename of pulse profile
%
% output parameters:
% Expt_out: updated Map with the experimental parameters
%
% April 2024 A. Kehl (akehl@gwdg.de)
%


function Expt_out=kehl_exp_pulse(expt_in,data,multipulses)

    % Check consistency
    grumble(expt_in,data,multipulses);
    Expt_out=expt_in;

    Expt_out("pulse")=data;
    Expt_out("3pulses")=multipulses;
end

function grumble(expt_in,data,multipulses)
if ~isa(expt_in,'containers.Map')
    error('expt_in must be a containers.Map object.');
end
if (~isnumeric(data))&&(~ischar(data))&&(~isstring(data))
    error('data must be numeric, a character string, or a string scalar.');
end
if (~islogical(multipulses))&&(~isnumeric(multipulses))
    error('multipulses must be logical, or numeric.');
end
end

