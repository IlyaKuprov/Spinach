% Static powder pulse-acquire experiment on a two-spin system 
% with a dipolar coupling.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function static_powder_dip()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={5.0,-2.0};
inter.coordinates={[0 0 0]; [0 3.9 0.1]};

% Algorithmic options
sys.disable={'trajlevel'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=1000;
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=12000;
parameters.npoints=128;
parameters.zerofill=512;
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

