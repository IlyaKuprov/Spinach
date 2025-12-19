% CPMG echo train in a powder.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function cpmg_echo_train()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5};
inter.zeeman.euler={[0 0 0],[0 0 0]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={ 50.0  50.0};
inter.r2_rates={150.0 150.0};
inter.equilibrium='zero';
inter.rlx_keep='secular';

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.grid='rep_2ang_200pts_sph';
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.pulse_op=operator(spin_system,'Lx','1H');
parameters.nloops=10;
parameters.timestep=1e-5;
parameters.npoints=100;
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@cpmg,parameters,'nmr');

% Plotting
x_axis=linspace(0,parameters.timestep*size(fid,2),size(fid,2));
kfigure(); plot(x_axis,real(fid)); kgrid; kxlabel('time, s');
kylabel('transverse magnetisation, a.u.'); xlim tight;

end

