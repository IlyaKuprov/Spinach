% PRESS excitation profile in three dimensions with tilted 
% gradient system. Change the frequency under Pulse Parame-
% ters to move the hot spot through the sample.
%
% Simulation time: hours, faster with a Tesla V100 GPU.
%
% a.j.allami@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function press_3d_example()

% Magnetic induction
sys.magnet=3.0;

% Spin systems
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={1.5,3.7};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=10;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Disable path tracing
sys.disable={'pt'};

% This is here essential
sys.enable={'polyadic'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.decouple={};
parameters.offset=0;
parameters.sweep=1000;
parameters.axis_units='ppm';
parameters.npoints=128;
parameters.invert_axis=1;
parameters.ss_grad_amp=[25e-3 25e-3 25e-3];
parameters.image_size=[101 103 105];    
parameters.sp_grad_amp=5e-3; % T/m
parameters.sp_grad_dur=5e-4;

% Pulse parameters
parameters.rf_phi={pi/2 pi/2 pi/2};
parameters.rf_frq_list={-30e3 +30e3 -50e3};
parameters.rf_amp_list={2*pi*5000 2*pi*5000 2*pi*5000};
parameters.rf_dur_list={0.5e-4 1.0e-4 1.0e-4};
parameters.max_rank={2 2 2};

% Sample geometry
parameters.dims=[0.30 0.25 0.27];
parameters.npts=[108 90 111];
parameters.deriv={'period',3};
parameters.grad_angles=[pi/3 pi/4 pi/5];

% Relaxation phantom
parameters.rlx_ph={zeros(parameters.npts)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts)};   
parameters.rho0_st={state(spin_system,'Lz','all')};
parameters.coil_ph={ones(parameters.npts)};
parameters.coil_st={state(spin_system,'L+','all')};

% Run voxel selection diagnostics
phan=imaging(spin_system,@press_voxel_3d,parameters);
volplot(phan,[-parameters.dims(1)/2 parameters.dims(1)/2 ...
              -parameters.dims(2)/2 parameters.dims(2)/2 ...
              -parameters.dims(3)/2 parameters.dims(3)/2]);
kxlabel('FOV1 / m'); kylabel('FOV2 / m'); kzlabel('FOV3 / m');

end

