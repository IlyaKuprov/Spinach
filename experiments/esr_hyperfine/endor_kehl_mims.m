%ENDOR_KEHL_MIMS Mims ENDOR pulse sequence for the SimDyn ENDOR context.
%
%   ENDOR_AMP=ENDOR_KEHL_MIMS(SPIN_SYSTEM,PARAMETERS,H,R,K) is called by
%   endor_kehl_context.m. H, R, and K are present for Spinach experiment
%   signature compatibility; this sequence uses the orientation-selected
%   effective Hamiltonian data prepared by the context.


function endor_amp=endor_kehl_mims(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);
if parameters.opt('Relax')==true
    endor_amp=kehl_mims_rlx(parameters.constants,parameters.spinOps,...
                                   parameters.spinSys,parameters.expt,parameters.opt,...
                                   parameters.paramsENDOR,parameters.epr);
else
    endor_amp=kehl_mims_calc(parameters.constants,parameters.spinOps,...
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

