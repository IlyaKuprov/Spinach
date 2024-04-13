% Powder NQR spectrum of a system with a single 127I nucleus.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function pure_nqr_iodine()

% System specification
sys.magnet=0;
sys.isotopes={'127I'};
inter.coupling.matrix{1,1}=eeqq2nqi(560e6,0.01,5/2,[0 0 0]);

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'damp'};
inter.damp_rate=1e5;
inter.rlx_keep='labframe';
inter.equilibrium='zero';
inter.temperature=298;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'127I'};
parameters.needs={'aniso_eq'};
parameters.sweep=5e8;
parameters.npoints=512;
parameters.grid='rep_2ang_200pts_sph';
parameters.coil=state(spin_system,'L+','127I');
parameters.pulse_op=(operator(spin_system,'L+','127I')+...
                     operator(spin_system,'L-','127I'))/2;
parameters.pulse_angle=pi/2;               
parameters.axis_units='MHz';

% Simulation
fid=powder(spin_system,@hp_acquire,parameters,'labframe');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spec=imag(fftshift(fft(fid)));

% Plotting
figure(); plot_1d(spin_system,spec,parameters);

end

