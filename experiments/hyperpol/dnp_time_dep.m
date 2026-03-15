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
% ilya.kuprov@weizmann.ac.il
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
if (~isnumeric(parameters.mw_pwr))||(~isreal(parameters.mw_pwr))||...
   (numel(parameters.mw_pwr)~=1)
    error('parameters.mw_pwr should be a real scalar.');
end
if ~isfield(parameters,'mw_off')
    error('microwave offset should be specified in parameters.mw_off variable.');
end
if (~isnumeric(parameters.mw_off))||(~isreal(parameters.mw_off))||...
   (numel(parameters.mw_off)~=1)
    error('parameters.mw_off should be a real scalar.');
end
if ~isfield(parameters,'rho0')
    error('thermal equilibrium state should be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('coil state should be specified in parameters.coil variable.');
end
if ~isfield(parameters,'mw_oper')
    error('microwave irradiation operator should be specified in parameters.mw_oper variable.');
end
if ~isfield(parameters,'ez_oper')
    error('electron Lz operator should be specified in parameters.ez_oper variable.');
end
if ~isfield(parameters,'nsteps')
    error('number of time steps should be specified in parameters.nsteps variable.');
end
if (~isnumeric(parameters.nsteps))||(~isreal(parameters.nsteps))||...
   (numel(parameters.nsteps)~=1)||(mod(parameters.nsteps,1)~=0)||...
   (parameters.nsteps<1)
    error('parameters.nsteps should be a positive real integer.');
end
if ~isfield(parameters,'dt')
    error('time step length should be specified in parameters.dt variable.');
end
if (~isnumeric(parameters.dt))||(~isreal(parameters.dt))||...
   (numel(parameters.dt)~=1)||(parameters.dt<=0)
    error('parameters.dt should be a positive real scalar.');
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

