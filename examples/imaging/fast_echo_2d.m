% Fast (in the experiment duration sense) spin echo 2D brain 
% imaging example.
%
% Simulation time: minutes, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function fast_echo_2d()

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
% sys.enable={'gpu'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=0.0;
parameters.image_size=[101 105];
parameters.ro_grad_amp=5.3e-3; % T/m
parameters.ro_grad_dur=2e-3;
parameters.pe_grad_amp=4.8e-3; % T/m
parameters.pe_grad_dur=1e-3;

% Relaxation superoperators
[R1Op,R2Op]=rlx_t1_t2(spin_system);

% Phantom library call
[R1Ph,R2Ph,PDPh,dims,npts]=phantoms('brain-medres');
R1Ph=R1Ph(:,:,50); R2Ph=R2Ph(:,:,50); PDPh=PDPh(:,:,50);

% Sample geometry
parameters.dims=dims([1 2]);
parameters.npts=npts([1 2]);
parameters.deriv={'period',3};

% Relaxation phantoms
parameters.rlx_ph={R1Ph,R2Ph};
parameters.rlx_op={R1Op,R2Op};

% Initial and detection state phantoms
parameters.rho0_ph={PDPh};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(parameters.npts)};
parameters.coil_st={state(spin_system,'L+','1H')};

% Run the simulation
mri=imaging(spin_system,@fse,parameters);

% Plotting
figure(); scale_figure([2.5 1.0]);
subplot(1,3,1); mri_2d_plot(mri,parameters,'image');    ktitle('recorded image');
subplot(1,3,2); mri_2d_plot(R1Ph,parameters,'phantom'); ktitle('$R_1$ phantom');
subplot(1,3,3); mri_2d_plot(R2Ph,parameters,'phantom'); ktitle('$R_2$ phantom');

end

