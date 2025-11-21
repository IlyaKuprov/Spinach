% Powder magic angle spinning spectrum of a single quadrupolar 
% deuterium nucleus using grid-free Fokker-Planck MAS formalism.
% Second order corrections to the rotating frame transformation
% are not applied.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function mas_powder_nqi_gridfree()

% System specification
sys.magnet=9.4;
sys.isotopes={'2H'};
inter.coupling.eigs={[-1e3 -2e3 3e3]};
inter.coupling.euler={[0.0 0.0 0.0]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.axis=[1 1 1];
parameters.max_rank=17;
parameters.rate=1e3;
parameters.sweep=2e4;
parameters.npoints=512;
parameters.zerofill=4096;
parameters.offset=0;
parameters.spins={'2H'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','2H','cheap');
parameters.coil=state(spin_system,'L+','2H','cheap');
parameters.verbose=1;

% Simulation
fid=gridfree(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

