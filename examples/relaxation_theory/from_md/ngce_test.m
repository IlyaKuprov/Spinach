% Test of the numerical integral route to the Redfield relaxation 
% superoperator against the analytical results for isotropic rota-
% tional diffusion.
%
% Calculation time: minutes.
%
% ilya.kuprov@weizmann.ac.il
% jpresteg@uga.edu

function ngce_test()

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
R_red=real(relaxation(spin_system));

% Get lab frame Hamiltonian components
[H0,Q]=hamiltonian(assume(spin_system,'labframe'));

% Get a random walk on a sphere
tau_c=inter.tau_c{1}; dt=tau_c/25;
eulers=rwalk(100000,tau_c,dt);

% Get Hamiltonian trajectory
H1=cell(size(eulers,1),1);
parfor n=1:size(eulers,1)
    H1{n}=orientation(Q,eulers(n,:));
end

% Compute the relaxation matrix using GCE
[R_gce,dR_gce]=ngce(spin_system,H0,H1,dt,tau_c,1e-3);

% Get the answers to the user
disp('Relaxation superoperator, analytical'); disp(full(R_red));
disp('Relaxation superoperator, numerical '); disp(full(R_gce));
disp('Relaxation superoperator, STD(mean) '); disp(full(dR_gce));

% Get relaxation rates
gce_rates_stdm=diag(dR_gce); 
gce_rates=diag(R_gce);
red_rates=diag(R_red);

% Draw a diagonal straight line
min_rate=min([gce_rates; red_rates]);
max_rate=max([gce_rates; red_rates]);
kfigure(); hold on;
plot([min_rate max_rate],...
     [min_rate max_rate],'b-');

% Plot 95% confidence bands for the rates
errorbar(gce_rates,red_rates,2*gce_rates_stdm,...
         'horizontal','LineStyle','none',...
         'Marker','o','Color','r'); hold on;
xlim padded; ylim padded; box on; kgrid;
kxlabel('numerical relaxation rates');
kylabel('analytical relaxation rates');
 
end

