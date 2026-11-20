% Stimulated-echo ESEEM with hard or rectangular electron pulses.
% The initial state is already transverse: the sequence is an outer
% delay, a 90-degree pulse, coherence selection, a variable storage
% delay, a second 90-degree pulse, and the same outer delay. Syntax:
%
%            fid=eseem_3p(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.rho0           - prepared state, column vector
%     parameters.coil           - detection state, column vector
%     parameters.pulse_op       - unit-amplitude RF commutator
%     parameters.pulse_duration - rectangular pulse duration, seconds;
%                                 zero selects instantaneous pulses
%     parameters.tau            - fixed outer delay, seconds
%     parameters.timestep       - storage delay increment, seconds
%     parameters.npoints        - number of storage delays, including zero
%
%     H, R, K - Hamiltonian, relaxation, and kinetics matrices from
%               the context, in a complete single-substance Liouville
%               basis (sphten-liouv or zeeman-liouv)
%
% Outputs:
%
%     fid - raw complex echo amplitudes as a column vector, using
%           the receiver phase supplied in parameters.coil
%
% Selection retains total electron order zero with total nuclear
% orders -1, 0, or +1. H, R, and K act during the free delays and
% during finite pulses, but not during instantaneous pulses. Supply the offset through
% the context, not the pulse sequence.
%
% Talos <talos@spindynamics.org>

function fid=eseem_3p(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Build the free and integrated pulse generators
L=H+1i*R+1i*K;
P=parameters.pulse_duration*L+(pi/2)*parameters.pulse_op;

% Prepare the selected storage state
rho=step(spin_system,L,parameters.rho0,parameters.tau);
rho=step(spin_system,P,rho,1);
rho=coherence(spin_system,rho,{{'electrons',0},{'nuclei',[-1 0 1]}});

% Pull the detector back through the last delay and pulse
coil=step(spin_system,-L',parameters.coil,parameters.tau);
coil=step(spin_system,-P',coil,1);

% Detect the echo as a function of storage delay
fid=evolution(spin_system,L,coil,rho,parameters.timestep,...
              parameters.npoints-1,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})||...
   (spin_system.bas.nsubst~=1)||...
   (spin_system.bas.offsets(end)~=prod(spin_system.comp.mults)^2)
    error('a complete single-substance Liouville basis is required.');
end
dim=spin_system.bas.offsets(end);
validateattributes(H,{'numeric'},{'size',[dim dim],'finite'});
validateattributes(R,{'numeric'},{'size',[dim dim],'finite'});
validateattributes(K,{'numeric'},{'size',[dim dim],'finite'});
validateattributes(parameters.rho0,{'numeric'},{'size',[dim 1],'finite'});
validateattributes(parameters.coil,{'numeric'},{'size',[dim 1],'finite'});
validateattributes(parameters.pulse_op,{'numeric'},{'size',[dim dim],'finite'});
validateattributes(parameters.pulse_duration,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(parameters.tau,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(parameters.timestep,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(parameters.npoints,{'numeric'},{'real','scalar','finite','integer','positive'});
end


