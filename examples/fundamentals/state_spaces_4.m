% Trajectory analysis for a MAS simulation of isotopically labelled 
% glycine powder, starting from L+ on protons.
%
% Calculation time: minutes on a Tesla A100 GPU.
%
% ilya.kuprov@weizmann.ac.il

function state_spaces_4()

% Spin system properties (PCM DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'),...
                     {{'H','1H'},{'C','13C'},{'N','15N'}},...
                     [31.5 182.1 264.5],[]);
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.longitudinals={'15N','13C'};

% Force Krylov propagation
sys.tols.krylov_tol=1000;
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=2000;
parameters.axis=[1 1 1];
parameters.max_rank=17;
parameters.sweep=1e5;
parameters.npoints=64;
parameters.offset=15000;
parameters.spins={'13C'};
parameters.grid='rep_2ang_200pts_sph';
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.verbose=1;

% Get the trajectory
traj=singlerot(spin_system,@traject,parameters,'nmr');

% Average over rotor phase
traj=fpl2rho(traj,2*parameters.max_rank+1);

% Trajectory analysis
figure(); trajan(spin_system,traj,'correlation_order');
xlim tight; ylim([1e-5 0.05]); set(gca,'YScale','log');

end

