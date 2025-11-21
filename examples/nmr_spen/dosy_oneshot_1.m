% Oneshot DOSY pulse sequence for a system of three coupled 
% spins with different relaxation rates.
%
% Timing: minutes on NVidia Tesla A100, much longer on CPU
%
% mariagrazia.concilio@sjtu.edu.cn

function dosy_oneshot_1()

% Magnetic field
sys.magnet=11.7428;    % Tesla

% Spin system
sys.isotopes={'1H','1H','1H'};
inter.zeeman.scalar={4.70 3.50 1.50};
inter.coupling.scalar{1,2}=15;
inter.coupling.scalar{2,3}=15;
inter.coupling.scalar{1,3}=10;
inter.coupling.scalar{3,3}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters 
inter.relaxation={'t1_t2'}; 
inter.rlx_keep='diagonal';  
inter.equilibrium='zero'; 
inter.r1_rates=num2cell(1./[0.1952 0.2100 0.2500]);
inter.r2_rates=num2cell(1./[0.1602 0.1802 0.1902]);

% Algorithmic options
sys.disable={'pt'};
sys.enable={'greedy'}; % 'gpu'

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Acquisition parameters
parameters.spins={'1H'};
parameters.sweep=5000;       % Hz
parameters.npoints=1024;  
parameters.zerofill=32768;
parameters.axis_units='ppm';
parameters.offset=2497.78;   % Hz

% Sample geometry
parameters.dims=0.015;       % m
parameters.npts=5000;
parameters.deriv={'period',7};

% Relaxation phantom
parameters.rlx_ph={ones(parameters.npts,1)};
parameters.rlx_op={relaxation(spin_system)};

% Initial and detection state phantoms
parameters.rho0_ph={ones(parameters.npts,1)};
parameters.rho0_st={state(spin_system,'Lz','1H','cheap')};
parameters.coil_ph={ones(parameters.npts,1)};
parameters.coil_st={state(spin_system,'L-','1H','cheap')};

% Diffusion
parameters.diff=18.55e-10;       % m^2/s

% Gradient parameters
parameters.g_amp=0.255;   % T/m 
parameters.kappa=0.2;         

% Timing parameters
parameters.g_dur=0.001;       % seconds
parameters.del=0.05;          % seconds
parameters.g_stab_del=0.0005; % seconds

% Call the sequence in the imaging context
fid=imaging(spin_system,@dosy_oneshot,parameters);

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',5}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,-real(spectrum),parameters);

end

