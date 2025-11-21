% Slice selection example using a one-dimensional sample and 
% a shaped slice selection pulse in the presence of diffusion
% and flow.
%
% Calculation time: seconds.
%
% Ahmed Allami
% Ilya Kuprov

function slice_select_1d_shaped()

% Isotopes
sys.isotopes={'1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={0.0};

% Relaxation model
inter.relaxation={'t1_t2'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.r1_rates={30};
inter.r2_rates={70};

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
parameters.axis_units='kHz';
parameters.invert_axis=1;
parameters.ss_grad_amp=30e-3;
parameters.ro_grad_amp=30e-3;

% Pulse parameters
pulse_nsteps=50; pulse_time=0.5e-4;
pulse_frq=+100e3; pulse_pwr=2*pi*20000;
parameters.rf_phi=pi/2;
parameters.rf_frq_list=pulse_frq*ones(1,pulse_nsteps);
parameters.rf_amp_list=pulse_pwr*pulse_shape('gaussian',pulse_nsteps);
parameters.rf_dur_list=(pulse_time/pulse_nsteps)*ones(1,pulse_nsteps);
parameters.max_rank=3;

% Sample geometry
parameters.dims=0.30;
parameters.npts=500;
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

% Run the pulse sequence in the imaging context
fid=imaging(spin_system,@slice_select_1d,parameters);

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'}});

% Run the Fourier transform
mri=real(fftshift(fft(ifftshift(fid))));

% Plotting
kfigure(); plot_1d(spin_system,mri,parameters);

end

