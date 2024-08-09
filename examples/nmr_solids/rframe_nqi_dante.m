% DANTE MAS spectrum of a single quadrupolar 14N nucleus using 1D 
% Fokker-Planck equation and a spherical grid. The calculation 
% accounts for the second-order quadrupolar shift and lineshape.
%
% Set to reproduce Figure 3d from
%
%           http://dx.doi.org/10.1016/j.jmr.2012.05.024
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk

function rframe_nqi_dante()

% System specification
sys.magnet=18.8; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.50,1,[0 0 0]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel'};
sys.enable={'prop_cache'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=62.5e3;
parameters.axis=[1 1 1];
parameters.max_rank=35;
parameters.grid='rep_2ang_200pts_sph';
parameters.sweep=2000000;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.offset=2200;
parameters.spins={'14N'};
parameters.rframes={{'14N',2}};
parameters.axis_units='Hz';
parameters.rho0=state(spin_system,'Lz','14N');
parameters.coil=state(spin_system,'L+','14N');
parameters.verbose=0;
parameters.pulse_dur=1.2e-6;
parameters.pulse_amp=88e3;
parameters.pulse_num=2;
parameters.n_periods=2;

% Simulation
fid=singlerot(spin_system,@dante,parameters,'labframe');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,abs(spectrum),parameters);

end

