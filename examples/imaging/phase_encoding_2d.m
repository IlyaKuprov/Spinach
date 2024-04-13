% Simple phase-encoded 2D imaging example.
%
% Calculation time: seconds.
%
% Ahmed Allami
% Ilya Kuprov

function phase_encoding_2d()

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
inter.r1_rates={30.0};
inter.r2_rates={70.0};

% Disable path tracing
sys.disable={'pt','krylov'};

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
parameters.ro_grad_amp=4.3e-3; % T/m
parameters.ro_grad_dur=2e-3;
parameters.pe_grad_amp=3.8e-3; % T/m
parameters.pe_grad_dur=1e-3;
parameters.t_echo=0.025;

% Sample geometry
parameters.dims=[0.30 0.25];
parameters.npts=[108 90];
parameters.deriv={'period',3};

% Relaxation phantom
[R1,R2]=rlx_t1_t2(spin_system);
load('../../etc/phantoms/letter_a.mat','R1Ph','R2Ph');
parameters.rlx_ph={R1Ph,R2Ph};
parameters.rlx_op={R1,R2};

% Initial and detection state phantoms
parameters.rho0_ph={ones(prod(parameters.npts,1))};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(prod(parameters.npts,1))};
parameters.coil_st={state(spin_system,'L+','1H')};

% Run the simulation
mri=imaging(spin_system,@phase_enc_2d,parameters);

% Plotting
loc=get(0,'defaultfigureposition'); figure('Position',[loc(1:2) 3*loc(3) loc(4)]);
subplot(1,3,1); mri_2d_plot(mri, parameters,'image');   ktitle('recorded image');
subplot(1,3,2); mri_2d_plot(R1Ph,parameters,'phantom'); ktitle('$R_1$ phantom');
subplot(1,3,3); mri_2d_plot(R2Ph,parameters,'phantom'); ktitle('$R_2$ phantom');

end

