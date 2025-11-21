% A gradient echo experiment in the presence
% of diffusion and flow.
%
% Calculation time: seconds.
%
% Ahmed Allami
% Ilya Kuprov

function gradient_echo_1d()

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
parameters.rlx_ph={};
parameters.rlx_op={};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H')};

% Diffusion and flow
parameters.u=ones(parameters.npts,1);
parameters.diff=1e-6;

% Run the simulation
echo=imaging(spin_system,@grad_echo,parameters);

% Plotting
grad_duration=parameters.g_step_dur*parameters.g_n_steps;
time_axis=linspace(-grad_duration,grad_duration,401);
kfigure(); plot(time_axis,real(echo)); axis('tight'); kgrid;
kxlabel('time, seconds'); kylabel('intensity, a.u.');
 
end

