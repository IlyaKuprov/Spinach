% A spin echo experiment under a gradient in the 
% presence of diffusion and flow.
%
% Calculation time: seconds.
%
% Ahmed Allami
% Ilya Kuprov

function spin_echo_1d()

% Isotopes
sys.isotopes={'1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={1.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.offset=0.0;
parameters.g_amp=5e-6;
parameters.g_step_dur=2e-4;
parameters.g_n_steps=200;

% Sample geometry
parameters.dims=0.30;
parameters.npts=100;
parameters.deriv={'period',3};

% Relaxation phantom
parameters.rlx_ph={zeros(parameters.npts,1)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H')};

% Diffusion and flow
parameters.u=ones(parameters.npts,1);
parameters.diff=1e-6;

% Run the simulation
echo=imaging(spin_system,@spin_echo,parameters);

% Plotting
time_axis=(0:400)*parameters.g_step_dur;
figure(); plot(time_axis',real(echo));
kylabel('echo intensity, a.u.');
kxlabel('time, seconds'); 
axis('tight'); kgrid;
 
end

