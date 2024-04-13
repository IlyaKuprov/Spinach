% Pulse-acquire NMR spectrum in a system with a hypothetical scalar
% coupling to a 235U nucleus. The spectral lines should be split
% accordingly.
%
% i.kuprov@soton.ac.uk

function high_spin_system_1()

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spin system
sys.isotopes={'1H','235U','1H','1H'};
inter.zeeman.scalar={-0.5  0.0  2.5  1.3};
inter.coupling.scalar{1,2}=100;
inter.coupling.scalar{3,4}=50;
inter.coupling.scalar{4,4}=0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Pulse sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=3500;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

