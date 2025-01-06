% DPFGSE water suppression example for a solution of GABA
% in water. Gradients and soft pulses are done explicitly.
%
% Simulation time: minutes, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% p.lally@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function dpfgse_signal_suppress()

% Isotopes
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H'};

% Magnet (Tesla)
sys.magnet=5.9;

% Chemical shifts (ppm)
inter.zeeman.scalar={3.00 3.00 1.88 1.88 2.28 2.28 4.80};

% J-couplings (Hz)
inter.coupling.scalar{1,3}=7.36;
inter.coupling.scalar{1,4}=7.36;
inter.coupling.scalar{2,3}=7.36;
inter.coupling.scalar{2,4}=7.36;
inter.coupling.scalar{3,5}=7.58;
inter.coupling.scalar{3,6}=7.58;
inter.coupling.scalar{4,5}=7.58;
inter.coupling.scalar{4,6}=7.58;
inter.coupling.scalar{7,7}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Disable path tracing
sys.disable={'pt','krylov'};

% This needs a GPU
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=800;
parameters.sweep=1200;
parameters.npoints=512;
parameters.zerofill=2048;
parameters.axis_units='Hz';
parameters.invert_axis=1;
parameters.g_amp=[1e-3 1.5e-3];
parameters.g_dur=1e-3;

% 180-degree shaped pulse on water
pulse_nsteps=10;
parameters.max_rank=2;
parameters.rf_phi=0;
parameters.rf_frq_list=1220*ones(1,pulse_nsteps);
parameters.rf_amp_list=2*pi*1700*pulse_shape('gaussian',pulse_nsteps);     
parameters.rf_dur_list=(2e-3/pulse_nsteps)*ones(1,pulse_nsteps);

% Sample geometry
parameters.dims=0.30;
parameters.npts=100;
parameters.deriv={'period',3};

% Relaxation phantoms and operators
parameters.rlx_ph={}; parameters.rlx_op={};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L+','1H')};

% No diffusion or flow
parameters.u=zeros(parameters.npts,1);
parameters.diff=0;

% Run simulation
fid=imaging(spin_system,@dpfgse_suppress,parameters);

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Run the Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

