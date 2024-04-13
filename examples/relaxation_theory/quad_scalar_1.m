% NMR spectrum of 17O enriched water inside a fullerene cage. A rather
% exotic combination of quadrupolar relaxation on the oxygen and H-O
% scalar coupling is driving proton relaxation in this case. 17O quad-
% rupolar parameters in gaseous (assumed to be) water are coming from
% http://dx.doi.org/10.1063/1.1672122 (Table III).
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function quad_scalar_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H','17O'};
inter.coupling.scalar=cell(3,3);
inter.coupling.scalar{1,2}=0;
inter.coupling.scalar{1,3}=80;
inter.coupling.scalar{2,3}=80;
inter.coupling.matrix{3,3}=eeqq2nqi(9.82e6,0.407,5/2,[0 0 0]);

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={1e-13};

% Basis specification
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
parameters.offset=0;
parameters.sweep=500;
parameters.npoints=512;
parameters.zerofill=2048;
parameters.axis_units='Hz';

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

