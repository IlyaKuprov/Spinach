% Powder magic angle spinning spectrum of a single anisotropically shielded
% proton spin using a Floquet theory based formalism.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function mas_powder_csa_floquet()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5};
inter.zeeman.euler={[0 0 0],[0 0 0]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=500;
parameters.axis=[1 1 1];
parameters.max_rank=17;
parameters.grid='leb_2ang_rank_17';
parameters.sweep=2e4;
parameters.npoints=512;
parameters.zerofill=4096;
parameters.offset=0;
parameters.spins={'1H'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.verbose=0;

% Simulation
fid=floquet(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

