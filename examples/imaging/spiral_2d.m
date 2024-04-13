% Spiral K-space imaging example in 2D.
%
% Calculation time: minutes.
%
% Ahmed Allami
% Ilya Kuprov

function spiral_2d()

%  Isotopes
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

% Algorithmic options
sys.disable={'pt'};

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
parameters.grad_amp=5e-2;
parameters.spiral_dur=5e-3;
parameters.spiral_npts=3000;
parameters.spiral_frq=5e4;
parameters.t_echo=0.050;

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
parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
parameters.coil_ph={ones(prod(parameters.npts,1))};
parameters.coil_st={state(spin_system,'L+','1H','cheap')};

% Diffusion and flow
parameters.u=zeros(parameters.npts);
parameters.v=zeros(parameters.npts);
parameters.diff=0;

% Run the simulation
mri=imaging(spin_system,@spiral,parameters);

% Translate gradients for plotting
parameters.ro_grad_dur=parameters.spiral_dur/sqrt(parameters.spiral_npts)/2;
parameters.pe_grad_dur=parameters.spiral_dur/sqrt(parameters.spiral_npts)/4;
parameters.ro_grad_amp=parameters.grad_amp;
parameters.pe_grad_amp=parameters.grad_amp;

% Plotting
loc=get(0,'defaultfigureposition'); figure('Position',[loc(1:2) 3*loc(3) loc(4)]);
subplot(1,3,1); mri_2d_plot(mri, parameters,'image');   title('recorded image');
subplot(1,3,2); mri_2d_plot(R1Ph,parameters,'phantom'); title('R1 phantom');
subplot(1,3,3); mri_2d_plot(R2Ph,parameters,'phantom'); title('R2 phantom');

end

