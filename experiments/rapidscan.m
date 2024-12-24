% Time-domain rapid field scan ESR experiment, Eatons style. Syntax:
%
%        [b_axis,spectrum]=rapidscan(spin_system,parameters)
%
% Parameters:
%
%    parameters.mw_pwr    - microwave power, rad/s
%
%    parameters.sweep     - magnetic field sweep extents,
%                           arouind the centre field specified
%                           in sys.magnet, a two-element vector
%                           in Tesla
%
%    parameters.nsteps    - number of steps in the magnetic 
%                           field sweep
%
%    parameters.timestep  - duration of each time step, seconds
%
% Outputs:
%
%    b_axis               - magnetic field axis, Tesla
%
%    spectrum             - L+ observable amplitude at each
%                           magnetic field
%
% Note: this experiment should be called directly without a context.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=rapidscan.m>

function [b_axis,spectrum]=rapidscan(spin_system,parameters)

% Check consistency
grumble(parameters);

% Get the microwave Hamiltonian
Ep=operator(spin_system,'L+','E');
H_mw=-parameters.mw_pwr*(Ep-Ep')/2i;

% Get the thermal equilibrium
H=hamiltonian(assume(spin_system,'labframe'),'left');
rho=equilibrium(spin_system,H);

% Relevant operators
Hz=hamiltonian(assume(spin_system,'labframe','zeeman'));
Hc=hamiltonian(assume(spin_system,'labframe','couplings'));
R=relaxation(spin_system);

% Put the drift Liouvillian together
L=Hz+Hc; L=(L+L')/2;

% Go into the microwave rotating frame
spin_system=assume(spin_system,'labframe');
C=carrier(spin_system,'E');
L=rotframe(spin_system,C,L,'E',1);

% Add microwave and relaxation terms
L=L+H_mw+1i*R;

% Get the detection state
coil=state(spin_system,'L+','E');

% Compute the waveform and the axis
waveform=linspace(parameters.sweep(1),...
                  parameters.sweep(2),parameters.nsteps);
b_axis=waveform+spin_system.inter.magnet;

% Normalize the Zeeman Hamiltonian
Hz=Hz/spin_system.inter.magnet;

% Preallocate the answer
spectrum=zeros([parameters.nsteps 1],'like',1i);
         
% Propagate the system
for n=1:parameters.nsteps
    spectrum(n)=coil'*rho;  
    rho=step(spin_system,L+waveform(n)*Hz,rho,parameters.timestep);
end

end

% Consistency enforcement
function grumble(parameters)
if ~isfield(parameters,'mw_pwr')
    error('microwave power must be specified in parameters.mw_pwr field.');
end
if (~isnumeric(parameters.mw_pwr))||(~isreal(parameters.mw_pwr))||...
   (~isscalar(parameters.mw_pwr))||(parameters.mw_pwr<0)
    error('parameters.mw_pwr must be a non-negative real scalar.');
end
if ~isfield(parameters,'sweep')
    error('sweep extents must be specified in parameters.sweep field.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))||...
   (numel(parameters.sweep)~=2)||(parameters.sweep(1)>=parameters.sweep(2))
    error('parameters.sweep must have two real elements in ascending order.');
end
if ~isfield(parameters,'nsteps')
    error('number of steps must be specified in parameters.nsteps field.');
end
if (~isnumeric(parameters.nsteps))||(~isreal(parameters.nsteps))||...
   (~isscalar(parameters.nsteps))||(mod(parameters.nsteps,1)~=0)||...
   (parameters.nsteps<1)
    error('parameters.nsteps must be a positive integer.');
end
if ~isfield(parameters,'timestep')
    error('time step must be specified in parameters.timestep field.');
end
if (~isnumeric(parameters.timestep))||(~isreal(parameters.timestep))||...
   (~isscalar(parameters.timestep))||(parameters.timestep<=0)
    error('parameters.mw_pwr must be a positive real scalar.');
end
end

% The smart way to keep people passive and obedient is to strictly
% limit the spectrum of acceptable opinion, but allow very lively
% debate within that spectrum.
%
% Noam Chomsky

