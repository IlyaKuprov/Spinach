% Slice selection in 3D followed by phase-encoded imaging 
% of the resulting slice. This simulation fills up a sys-
% tem with eight H200 GPUs and 4 TB of RAM.
%
% Simulation time: you hope and pray this even starts;
%                  if it does, then hours.
%
% ilya.kuprov@weizmann.ac.il

function ph_enc_3d_highres()

% Isotopes
sys.isotopes={'1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={0.0};

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={1.0};
inter.r2_rates={1.0};
inter.equilibrium='zero';
inter.rlx_keep='diagonal';

% Disable path tracing
sys.disable={'pt'};

% This needs a GPU
sys.enable={'greedy','polyadic','gpu'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Gat phantom from library
[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-highres');

% Sample settings
parameters.dims=dims;
parameters.npts=npts;
parameters.deriv={'period',3};
parameters.spins={'1H'};
parameters.offset=0.0;

% Relaxation operators and phantoms
[R1,R2]=rlx_t1_t2(spin_system);
parameters.rlx_op={R1,R2};
parameters.rlx_ph={R1_Ph,R2_Ph};

% Initial and detection state phantoms
parameters.rho0_ph={PD_Ph};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(prod(parameters.npts,1))};
parameters.coil_st={state(spin_system,'L+','1H')};

% Diffusion and flow
parameters.u=zeros(parameters.npts);
parameters.v=zeros(parameters.npts);
parameters.w=zeros(parameters.npts);
parameters.diff=0;

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

% Imaging parameters
parameters.image_size=[129 129];
parameters.ss_grad_amp=32.0e-3;  % T/m
parameters.ss_grad_dur=1.0e-4;
parameters.pe_grad_amp=32.0e-3;  % T/m
parameters.pe_grad_dur=2.0e-4;
parameters.ro_grad_amp=32.0e-3;  % T/m
parameters.ro_grad_dur=3.0e-4;
parameters.t_echo=20e-3;
parameters.grad_angles=[pi/3 pi/4 pi/5];

% Draw phantom at high contrast
figure(); dims=zeros(1,6); 
dims([1 3 5])=-parameters.dims/2;
dims([2,4,6])=+parameters.dims/2;
volplot(R1_Ph,dims);
ktitle('$R_1$ map in 3D'); drawnow();

% Slice image using phase encoding sequence
fid=imaging(spin_system,@phase_enc_3d,parameters);

% Plot k-space representation of the slice
figure(); scale_figure([2.0 1.0]); subplot(1,2,1);
mri_2d_plot(fid.^(1/4),parameters,'k-space');
ktitle('$k$-space representation');

% Apodisation
fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}});

% Fourier transform
mri=real(fftshift(fft2(ifftshift(fid))));

% Plot real space representation of the slice
subplot(1,2,2); mri_2d_plot(mri,parameters,'image');
ktitle('real space representation');
                 
end

