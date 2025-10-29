% 2D diffusion weighted image with an arbitrary geometric
% pattern serving as diffusion coefficient distribution.
% 
% Simulation time: minutes, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function diffusion_weighted_2d()

% Isotopes
sys.isotopes={'1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={0.0};

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
parameters.diff_g_amp=[1e-3 1e-3];    % T/m
parameters.ro_grad_amp=4.3e-3; % T/m
parameters.ro_grad_dur=2e-3;
parameters.pe_grad_amp=3.8e-3; % T/m
parameters.pe_grad_dur=1e-3;
parameters.t_echo=1e-2;

% Sample geometry
parameters.dims=[0.30 0.25];
parameters.npts=[90 108];
parameters.deriv={'period',3};

% Relaxation phantoms and operators
parameters.rlx_ph={}; parameters.rlx_op={};

% Initial and detection state phantoms
parameters.rho0_ph={ones(prod(parameters.npts,1))};
parameters.rho0_st={state(spin_system,'Lz','1H')};
parameters.coil_ph={ones(prod(parameters.npts,1))};
parameters.coil_st={state(spin_system,'L+','1H')};

% Diffusion and flow
parameters.u=zeros(parameters.npts);
parameters.v=zeros(parameters.npts);

% 2D diffusion tensor field
load('../../etc/phantoms/pattern.mat','pattern');
parameters.dxx=1e-3*pattern;
parameters.dxy=zeros(parameters.npts);
parameters.dyx=zeros(parameters.npts);
parameters.dyy=1e-3*pattern;
 
% Run the simulation
mri=imaging(spin_system,@phase_enc_2d,parameters);

% Plotting
loc=get(0,'defaultfigureposition'); 
figure('Position',[loc(1:2) 2*loc(3) loc(4)]);
subplot(1,2,1); mri_2d_plot(mri,parameters,'image'); 
ktitle('diff.-weighted image');
subplot(1,2,2); mri_2d_plot(pattern,parameters,'phantom'); 
ktitle('diff. coeff. distribution');

end

