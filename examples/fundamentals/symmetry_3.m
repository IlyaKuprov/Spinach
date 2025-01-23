% 1H NMR spectrum of valine. Uses the fully symmetric 
% irreducible representation of the complete S6 group.
%
% ilya.kuprov@weizmann.ac.il

function symmetry_3()

% System specificaton
sys.magnet=11.7;
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'};
inter.zeeman.scalar={3.5950 2.2580 1.0270 1.0270 1.0270 0.9760 0.9760 0.9760};
inter.coupling.scalar{1,2}=4.34;
inter.coupling.scalar{2,3}=7.00;
inter.coupling.scalar{2,4}=7.00;
inter.coupling.scalar{2,5}=7.00;
inter.coupling.scalar{2,6}=7.00;
inter.coupling.scalar{2,7}=7.00;
inter.coupling.scalar{2,8}=7.00;
inter.coupling.scalar{8,8}=0.00;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.sym_spins={[3 4 5 6 7 8]};
bas.sym_group={'S6'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=1000;
parameters.sweep=2500;
parameters.npoints=8192;
parameters.zerofill=65536;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',5}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

