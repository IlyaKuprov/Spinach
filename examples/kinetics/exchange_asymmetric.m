% Two-spin asymmetric chemical exchange pattern.
%
% Calculation time: seconds.
%
% i.kuprov@soton.ac.uk

function exchange_asymmetric()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,3};
inter.chem.parts={1,2};
inter.chem.rates=[-5e2  2e3; 
                   5e2 -2e3];
inter.chem.concs=[2e3 5e2];

% Basis specification
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','chem');
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

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

