% Cross-polarisation experiment in the rotating frame. Applies a
% soft pi/2 pulse using the specified operators, then evolves the
% system with the specified spin-lock terms added to the Liovilli-
% an. The contact curve is returned. Syntax:
%
%   contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.spins      - working spins, a cell array of
%                             strings with high-gamma spins
%                             first and low-gamma spins last,
%                             for example {'1H','13C'}
%
%     parameters.hi_pwr     - nutation frequency of the exci-
%                             tation pulse on the high-gamma
%                             spins, Hz
%
%     parameters.cp_pwr     - nutation frequencies on the two
%                             channels during the CP contact
%                             time, a two-element vector, Hz
%
%     parameters.timestep   - time step of the CP contact ti-
%                             me, seconds
%
%     parameters.nsteps     - number of time steps to take
%                             during the CP contact time
%
%     parameters.rho0       - initial state, the state of the
%                             low-gamma spins will be wiped
%
%     parameters.coil       - detection state vector
%
%     H - Hamiltonian matrix, received from context function
%
%     R - relaxation superoperator, received from context function
%
%     K - kinetics superoperator, received from context function
%
% Output:
%
%     contact_curve - contact curve detected on the coil
%                     state specified in parameters.coil
%
% marina.carravetta@soton.ac.uk
% p.t.williamson@soton.ac.uk
% guinevere.mathies@uni-konstanz.de
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=cp_contact_soft.m>

function contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)

% Check consistency
grumble(parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Wipe the state of 13C
[~,rho]=decouple(spin_system,[],parameters.rho0,...
                                parameters.spins(2));

% Build and project 1H and 13C control operators
Hx=operator(spin_system,'Lx',parameters.spins{1});
Hy=operator(spin_system,'Ly',parameters.spins{1});
Cx=operator(spin_system,'Lx',parameters.spins{2});
Hx=kron(speye(parameters.spc_dim),Hx);
Hy=kron(speye(parameters.spc_dim),Hy);
Cx=kron(speye(parameters.spc_dim),Cx);

% Apply the 90-degree pulse on 1H along +X
rho=step(spin_system,L+2*pi*parameters.hi_pwr*Hx,...
                     rho,1/(4*parameters.hi_pwr));

% Run the CP contact time evolution: irradiation 
% of 1H along -Y, and of 13C along +X 
contact_curve=evolution(spin_system,L-2*pi*parameters.cp_pwr(1)*Hy...
                                     +2*pi*parameters.cp_pwr(2)*Cx,...
                        parameters.coil,rho,parameters.timestep,...
                        parameters.nsteps,'observable');

end

% Consistency enforcement
function grumble(parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if ~isfield(parameters,'hi_pwr')||(parameters.hi_pwr<=0)
    error('high RF amplitude must be specified in parameters.hi_pwr variable.');
end
if ~isfield(parameters,'cp_pwr')
    error('CP RF amplitudes must be specified in parameters.cp_pwr variable.');
end
if ~isfield(parameters,'rho0')
    error('initial state must be specified in parameters.rho0 variable.');
end
if ~isfield(parameters,'coil')
    error('detection state must be specified in parameters.coil variable.');
end
if ~isfield(parameters,'nsteps')
    error('number of time steps must be specified in parameters.nsteps variable.');
end
if (~isnumeric(parameters.nsteps))||(numel(parameters.nsteps)~=1)||...
   (~isreal(parameters.nsteps))||(parameters.nsteps<1)||(mod(parameters.nsteps,1)~=0)
    error('parameters.nsteps must be a positive integer.');
end
if ~isfield(parameters,'timestep')
    error('time step must be specified in parameters.timestep variable.');
end
if (~isnumeric(parameters.timestep))||(~isscalar(parameters.timestep))||...
   (~isreal(parameters.timestep))||(parameters.timestep<0)
    error('parameters.timestep must be a positive real scalar.');
end
end

% Those who lack the courage will always find a 
% philosophy to justify it.
%
% Albert Camus

