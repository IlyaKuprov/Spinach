% Time-domain spin dynamics under microwave irradiation. Syntax:
%
%        answer=dnp_time_dep(spin_system,parameters,H,R,K)
%  
% Parameters:
%
%    parameters.mw_pwr       -  microwave power, rad/s
%
%    parameters.mw_off       -  microwave frequency offset
%                               from free electron,rad/s
%
%    parameters.rho0         -  thermal equilibrium state
%
%    parameters.coil         -  coil state vector or a hori-
%                               zontal stack thereof
%
%    parameters.mw_oper      -  microwave irradiation operator
%
%    parameters.ez_oper      -  Lz operator on the electrons
%
%    parameters.dt           -  time step, seconds
%
%    parameters.nsteps       -  number of time steps
%
%    H - Hamiltonian matrix, received from context function
%
%    R - relaxation superoperator, received from context function
%
%    K - kinetics superoperator, received from context function
%
% Outputs:
%
%    answer  - a matrix of projections of the trajectory on 
%              each of the coils provided at each time step
%
% Note: the relaxation superoperator must be thermalised for this 
%       type of calculation.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dnp_time_dep.m>

function answer=dnp_time_dep(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters)

% Add microwave terms to the Hamiltonian
H=H+parameters.mw_pwr*parameters.mw_oper;

% Add microwave offset to the Hamiltonian
H=H-parameters.mw_off*parameters.ez_oper;

% Run the time evolution
answer=evolution(spin_system,H+1i*R+1i*K,parameters.coil,parameters.rho0,...
                 parameters.dt,parameters.nsteps,'multichannel');

end

% Consistency enforcement
function grumble(spin_system,parameters)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available for sphten-liouv and zeeman-liouv formalisms.');
end
if ~isfield(parameters,'mw_pwr')
    error('microwave power should be specified in parameters.mw_pwr variable.');
end
if numel(parameters.mw_pwr)~=1
    error('parameters.mw_pwr should have exactly one element.');
end
if ~isfield(parameters,'nsteps')
    error('number of time steps should be specified in parameters.nsteps variable.');
end
if numel(parameters.nsteps)~=1
    error('parameters.nsteps should have exactly one element.');
end
if ~isfield(parameters,'dt')
    error('time step length should be specified in parameters.dt variable.');
end
if numel(parameters.dt)~=1
    error('parameters.dt should have exactly one element.');
end
end

% A pack of wolves chased and cornered a shepherd dog who, faced with 
% the alternative of certain death, agreed to assist them in stealing
% sheep. The dog was accepted into the pack and granted all the atten-
% dant privileges. Three years on, during a particularly tough winter,
% the wolves had no choice but to eat the dog. The skin and the bones
% were given a proper burial; the wolves felt guilty. They spent some
% time deciding the epitaph: "here lies a friend" was not right since
% they ate the dog; neither was he a foe. They wrote: "colleague".
%
% Russian internet folklore, 
% translated by IK

