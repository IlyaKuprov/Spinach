% A Primas-style stochastic NMR experiment on strychnine. The calculation
% requires an NVidia Titan V GPU at a minimum.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function snmr_strychnine()

% Read the spin system properties 
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=5.9;

% Algorithmic options
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='IME';
inter.temperature=298;
inter.rlx_keep='kite';
inter.tau_c={200e-12};

% Use GPU arithmetic
% sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get the Hamiltonian
H0=hamiltonian(assume(spin_system,'nmr'));

% Get the relaxation superoperator
R=relaxation(spin_system);

% Stochastic process parameters
dt=1e-5;       % seconds 
omega_max=1e2; % rad/s
nsteps=10000;

% Control operators
Hx=operator(spin_system,'Lx','1H');
Hy=operator(spin_system,'Ly','1H');

% Uniformly distributed noise
Cx=omega_max*(2*rand(nsteps,1)-1);
Cy=omega_max*(2*rand(nsteps,1)-1);

% Isotropic thermal equilibrium
rho_eq=equilibrium(spin_system);

% GPU uploads
H0=gpuArray(H0); R=gpuArray(R);
Hx=gpuArray(Hx); Hy=gpuArray(Hy);

% Trajectory calculation
report(spin_system,'Computing trajectory...');
traj=gpuArray.zeros(numel(rho_eq),nsteps+1); 
traj(:,1)=gpuArray(rho_eq); tic;
for n=1:nsteps

    % Evolution generator
    G=H0+1i*R+Cx(n)*Hx+Cy(n)*Hy;
    
    % Time step
    traj(:,n+1)=step(spin_system,G,traj(:,n),dt);

end

% Performance report
disp(['Steps per second: ' num2str(nsteps/toc)]);

% GPU download
traj=gather(traj);

% Observables
coil=state(spin_system,'L+','1H'); fid=coil'*traj;

% Plotting - control sequence
kfigure(); subplot(1,2,1);
time_axis=linspace(0,nsteps*dt,nsteps);
plot(time_axis,[Cx Cy]');
axis tight; kgrid; title('controls');
legend({'Cx','Cy'},'Location','NorthEast');
kxlabel('time, seconds'); 
kylabel('nutation frequency, rad/s');

% Plotting - trajectory
subplot(1,2,2); scale_figure([1.5 1]);
time_axis=linspace(0,nsteps*dt,nsteps+1);
plot(time_axis,[real(fid); imag(fid)]);
axis tight; kgrid; title('trajectory');
legend({'Hx','Hy'},'Location','NorthEast');
kxlabel('time, seconds');
kylabel('absolute polarisation');

end

