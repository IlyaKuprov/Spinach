% Transverse magnetisation dynamics in a pulse-acquire experiment on
% a fatty acid. This example looks at how the magnetisation drifts
% around the state space under the influence of strong J-coupling in 
% the absence of relaxation.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk

function state_spaces_3()

% Read spin system properties 
[sys,inter]=fatty_acid(15);

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Algorithmic options
sys.enable={'greedy','prop_cache'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'Lx','1H');
parameters.coil=state(spin_system,'L+','1H');

% Assumptions
spin_system=assume(spin_system,'nmr');

% Operators
H=hamiltonian(spin_system);
Lx=operator(spin_system,'Lx','1H');

% Trajectory generation
traj=evolution(spin_system,H,[],parameters.rho0,4e-5,50,'trajectory');

% Number of repetition
N_rep=8;

% CPMG loop
for n=1:N_rep
    
    % Pulse
    traj(:,end)=step(spin_system,Lx,traj(:,end),pi);
    
    % Trajectory generation
    traj=[traj evolution(spin_system,H,[],traj(:,end),4e-5,100,'trajectory')]; %#ok<AGROW>
    
end

% Trajectory analysis
figure(); trajan(spin_system,traj,'correlation_order');

end

