% Test of the numerical integral route to the Redfield relaxation 
% superoperator against the analytical results for isotropic rota-
% tional diffusion.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% jpresteg@uga.edu

function ngce_test()

rng('shuffle');

% System specification
sys.magnet=9.4;
sys.isotopes={'1H','13C'};
inter.coordinates={[0.0 0.0 0.0]
                   [0.0 0.0 1.02]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1.0e-10};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get Redfield relaxation matrix
R_red=relaxation(spin_system);

% Get Hamiltonian components
[H0,Q]=hamiltonian(assume(spin_system,'labframe'));

% Get a random walk on a sphere
tau_c=inter.tau_c{1}; dt=tau_c/25;
eulers=rwalk(100000,tau_c,dt);

% Get Hamiltonian trajectories
H1=cell(size(eulers,1),1);
parfor n=1:size(eulers,1)
    H1{n}=orientation(Q,eulers(n,:));
end

% Get the GCE relaxation matrix
R_gce=ngce(spin_system,H0,H1,dt,tau_c,1e-3);

% Get the answers to the user
disp('Relaxation superoperator, numerical'); disp(full(R_gce));
disp('Relaxation superoperator, analytical'); disp(full(R_red));

% Do diagnostic plotting
figure();
gce_rates=diag(R_gce);
red_rates=diag(R_red);
min_rate=min([R_gce; R_red]);
max_rate=max([R_gce; R_red]);
plot(gce_rates,red_rates,'ro'); hold on;
plot([min_rate max_rate],[min_rate max_rate],'b-');
axis tight; box on; kgrid;
kxlabel('numerical relaxation rates');
kylabel('analytical relaxation rates');

% Get an accuracy indicator
disp(['RMSD numerical vs analytical: ' ... 
      num2str(sqrt(mean((gce_rates-red_rates).^2))) ' Hz']);
  
end

