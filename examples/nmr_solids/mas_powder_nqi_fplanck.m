% Powder magic angle spinning spectrum of a single quadrupolar
% deuterium nucleus using Fokker-Planck theory. Perturbative cor-
% rections to the rotationg frame transformation are not applied.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function mas_powder_nqi_fplanck()

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
parameters.rate=1000;
parameters.axis=[1 1 1];
parameters.max_rank=17;
parameters.grid='leb_2ang_rank_17';
parameters.sweep=2e4;
parameters.npoints=512;
parameters.zerofill=4096;
parameters.offset=0;
parameters.spins={'2H'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','2H');
parameters.coil=state(spin_system,'L+','2H');
parameters.verbose=0;

% Simulation
fid=singlerot(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

