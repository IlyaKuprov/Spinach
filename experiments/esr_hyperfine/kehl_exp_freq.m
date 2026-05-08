% adds experimental parameters for freq domain EPR calculations
%
% input parameters:
% expt_in: the Map with the experimental parameters
% freq_min: minimal sweep freq in GHz
% freq_range: range of freq sweep in GHz
% freq_steps: steps of freq sweep
%
% output parameters:
% Expt_out: updated Map with the experimental parameters
%
% February 2024 A. Kehl (akehl@gwdg.de)
%

function Expt_out=kehl_exp_freq(expt_in,freq_min,freq_range,freq_steps)

    % Check consistency
    grumble(expt_in,freq_min,freq_range,freq_steps);
    Expt_out=expt_in;

    Expt_out("FreqMin")=freq_min*1e9;
    Expt_out("FreqRange")=freq_range*1e9;
    Expt_out("FreqSteps")=freq_steps*1e9;

end

function grumble(expt_in,freq_min,freq_range,freq_steps)
if ~isa(expt_in,'containers.Map')
    error('expt_in must be a containers.Map object.');
end
if ~isnumeric(freq_min)
    error('freq_min must be numeric.');
end
if ~isnumeric(freq_range)
    error('freq_range must be numeric.');
end
if ~isnumeric(freq_steps)
    error('freq_steps must be numeric.');
end
end

