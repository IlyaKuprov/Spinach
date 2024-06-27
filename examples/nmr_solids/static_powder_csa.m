% Two-spin static CSA powder pattern.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function static_powder_csa()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5};
inter.zeeman.euler={[0 0 0],[0 0 0]};

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
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=15000;
parameters.npoints=256;
parameters.zerofill=512;
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

