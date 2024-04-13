% Two-spin asymmetric magnetization flux problem.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk

function flux_asymmetric()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,3};
inter.chem.flux_rate=zeros(2);
inter.chem.flux_rate(1,2)=5e2; % from 1 to 2
inter.chem.flux_rate(2,1)=2e3; % from 2 to 1
inter.chem.flux_type='intermolecular';

% Basis specification
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=2e3*state(spin_system,{'L+'},{1})+...  % Initial concentrations
                5e2*state(spin_system,{'L+'},{2});     % are specified here
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=900;
parameters.sweep=5000;
parameters.npoints=512;
parameters.zerofill=1024;
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

