% 1D imaging experiment with a hard pulse in 
% the presence of diffusion and flow.
%
% Calculation time: seconds.
%
% Ahmed Allami
% Ilya Kuprov

function gradient_image_1d()

% Isotopes
sys.isotopes={'1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={1.0};

% Relaxation model
inter.relaxation={'t1_t2'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.r1_rates={0.5};
inter.r2_rates={2.0};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=0;
parameters.sweep=500000;
parameters.npoints=128;
parameters.zerofill=128;
parameters.axis_units='kHz';
parameters.invert_axis=1;
parameters.ro_grad_amp=30e-3;

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
parameters.u=1e-2*ones(parameters.npts,1);
parameters.diff=5e-6;

% Run the simulation
fid=imaging(spin_system,@basic_1d_hard,parameters);

% Apodisation
fid=apodization(fid,'sqsinbell-1d');

% Run the Fourier transform
mri=real(fftshift(fft(ifftshift(fid))));

% Plotting
figure(); plot_1d(spin_system,mri,parameters);

end

