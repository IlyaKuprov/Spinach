% Slice selection in 3D followed by three-dimensional echo
% planar imaging sequence.
%
% Simulation time: minutes, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function echo_planar_3d()

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
inter.r1_rates={1};
inter.r2_rates={1};

% Disable path tracing
sys.disable={'pt'};

% This needs a GPU
sys.enable={'greedy'}; % 'gpu'

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Pulse phase
parameters.rf_phi=pi/2;

% Number of steps in the pulse
pulse_nsteps=50;

% Overall pulse duration
pulse_time=2.0e-4;

% Slice selection pulse frequency table
parameters.rf_frq_list=-5e3*ones(1,pulse_nsteps);

% Slice selection pulse amplitude table
parameters.rf_amp_list=2*pi*7500*pulse_shape('gaussian',pulse_nsteps);

% Slice selection pulse duration table
parameters.rf_dur_list=(pulse_time/pulse_nsteps)*ones(1,pulse_nsteps);

% Sequence parameters
parameters.image_size=[201 201];
parameters.ss_grad_amp=32.0e-3;  % T/m
parameters.ss_grad_dur=1.0e-3;
parameters.ro_grad_amp=5.3e-3; % T/m
parameters.pe_grad_amp=4.8e-3; % T/m
parameters.ro_grad_dur=4e-3;
parameters.pe_grad_dur=4e-3;
parameters.t_echo=20e-3;
parameters.grad_angles=[pi/3 pi/4 pi/5];

% Phantom library call
[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-medres');

% Sample settings
parameters.dims=dims;
parameters.npts=npts;
parameters.deriv={'period',3};
parameters.spins={'1H'};
parameters.offset=0.0;

% Relaxation phantom
[R1,R2]=rlx_t1_t2(spin_system);
parameters.rlx_op={R1,R2};
parameters.rlx_ph={R1_Ph,R2_Ph};

% Initial and detection state phantoms
parameters.rho0_ph={PD_Ph};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(prod(parameters.npts,1))};
parameters.coil_st={state(spin_system,'L+','1H')};

% No diffusion or flow
parameters.u=zeros(parameters.npts);
parameters.v=zeros(parameters.npts);
parameters.diff=0;

% Draw the phantom
figure(); dims=zeros(1,6); 
dims([1 3 5])=-parameters.dims/2;
dims([2 4 6])=+parameters.dims/2;
volplot(R1_Ph,dims);
ktitle('$R_1$ map in 3D'); drawnow();

% Run the simulation
fid=imaging(spin_system,@epi_3d,parameters);

% Plot k-space representation of the slice,
% the .^(1/4) improves fringe visibility
figure(); scale_figure([2.0 1.0]); subplot(1,2,1);
mri_2d_plot(fid.^(1/4),parameters,'k-space');
ktitle('$k$-space representation');

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}});

% Fourier transform
mri=real(fftshift(fft2(ifftshift(fid))));

% For FOV calculation, G{1} is effectively halved
parameters.pe_grad_amp=parameters.pe_grad_amp/2;

% Plot real space representation of the slice
subplot(1,2,2); mri_2d_plot(mri,parameters,'image');
ktitle('real space representation');

end

