%ENDOR_KEHL_TIME Time-domain Mims ENDOR sequence for endor_kehl_context.m.


function endor_amp=endor_kehl_time(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.opt('Relax')==true
    endor_amp=kehl_time_rlx(parameters.constants,parameters.spinOps,...
                                         parameters.spinSys,parameters.expt,parameters.opt,...
                                         parameters.paramsENDOR,parameters.epr);
else
    endor_amp=kehl_time_calc(parameters.constants,parameters.spinOps,...
                                   parameters.spinSys,parameters.expt,parameters.opt,...
                                   parameters.paramsENDOR,parameters.epr);
end
end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
    error('spin_system must be a Spinach spin system structure.');
end
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
if (~isempty(H))&&(~isnumeric(H))
    error('H must be empty or numeric.');
end
if (~isempty(R))&&(~isnumeric(R))
    error('R must be empty or numeric.');
end
if (~isempty(K))&&(~isnumeric(K))
    error('K must be empty or numeric.');
end
end

