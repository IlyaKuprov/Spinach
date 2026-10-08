% Powder magic angle spinning spectrum (rotor-synchronized detection)
% of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation
% and a spherical grid. The calculation accounts for the second-order
% quadrupolar shift and lineshape by applying numerical second order
% corrections to the rotating frame transformation.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function rframe_nqi_fplanck()

% System specification
sys.magnet=14.1; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(3.06e6,0.40,1,[0 0 0]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'trajlevel','krylov'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=50000;
parameters.axis=[1 1 1];
parameters.max_rank=85;
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=50000;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.offset=18000;
parameters.spins={'14N'};
parameters.rframes={{'14N',2}};
parameters.axis_units='Hz';
parameters.rho0=state(spin_system,'L+','14N');
parameters.coil=state(spin_system,'L+','14N');
parameters.verbose=0;

% Simulation
fid=singlerot(spin_system,@acquire,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

