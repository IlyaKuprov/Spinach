% A garden variety strongly coupled two-spin system.
%
% i.kuprov@soton.ac.uk

function strong_coupling()

% Isotopes
sys.isotopes={'1H','1H'};

% Magnetic induction
sys.magnet=5.9;

% Chemical shifts
inter.zeeman.scalar={1.0 1.5};

% Scalar couplings
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=7.0; 

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=300;
parameters.sweep=300;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',10);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

