% Correlation order dynamics in a pulse-acquire experiment on
% strychnine. Set to reproduce Figure 4 from our state space
% restriction accuracy analysis paper:
%
%             http://dx.doi.org/10.1063/1.3624564
%
% Run time: hours (much faster on a Tesla A100 GPU)
%
% ilya.kuprov@weizmann.ac.il

function state_spaces_1()

% Read spin system properties 
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.level=7; bas.space_level=1;
bas.connectivity='scalar_couplings';
bas.projections=1;

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Algorithmic options
sys.disable={'trajlevel'};
sys.enable={'greedy'}; % 'gpu'

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='kite';
inter.tau_c={200e-12};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Initial condition
rho=state(spin_system,'L+','1H');

% Assumptions
spin_system=assume(spin_system,'nmr');

% Liouvillian
L=hamiltonian(spin_system)+1i*relaxation(spin_system);

% Trajectory generation
traj=evolution(spin_system,L,[],rho,1e-3,1000,'trajectory');

% Trajectory analysis
figure(); trajan(spin_system,traj,'correlation_order'); 
axis([0 1000 1e-7 10]); set(gca,'YScale','log');

end

