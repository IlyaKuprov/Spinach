% Double angle spinning spectrum of N-acetylvaline 14N nucleus
% using 1D Fokker-Planck equation and a spherical grid. The cal-
% culation includes the second-order quadrupolar shift and the
% third-order lineshape. Time-domain detection.
%
% Note: slower spinning rates and larger NQIs require larger 
%       ranks and spherical grids. At the moment the spinning
%       frequencies are set artificially too high to reduce
%       the simulation time in this example.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function dor_powder_nav_fplanck_time()

% System specification
sys.magnet=14.1; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(3.21e6,0.27,1,[0 0 0]);

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=2e3;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate_outer=1e6;
parameters.rate_inner=5e6;
parameters.rank_outer=7;
parameters.rank_inner=4;
parameters.axis_outer=[sqrt(2/3) 0 sqrt(1/3)];                     % 54.74 degrees
parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))]; % 30.56 degrees
parameters.grid='rep_2ang_100pts_oct';
parameters.sweep=1e5;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.spins={'14N'};
parameters.axis_units='kHz';
parameters.rho0=state(spin_system,'L+','14N','cheap');
parameters.coil=state(spin_system,'L+','14N','cheap');
parameters.rframes={{'14N',3}};
parameters.verbose=0;

% Simulation
fid=doublerot(spin_system,@acquire,parameters,'labframe');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

