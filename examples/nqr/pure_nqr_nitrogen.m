% Powder NQR spectrum of a system with a single 14N nucleus.
%
% Calculation time: seconds
%
% lewis.robertson@csiro.au
% ilya.kuprov@weizmann.ac.il

function pure_nqr_nitrogen()

% System specification
sys.magnet=0;
sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);

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
parameters.spins={'14N'};
parameters.needs={'aniso_eq'};
parameters.sweep=5e6;
parameters.npoints=512;
parameters.grid='rep_2ang_200pts_sph';
parameters.coil=state(spin_system,'L+','14N');
parameters.pulse_op=operator(spin_system,'Lx','14N');
parameters.pulse_angle=pi/2;               
parameters.axis_units='MHz';

% Simulation
fid=powder(spin_system,@hp_acquire,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spec=imag(fftshift(fft(fid)));

% Plotting
kfigure(); plot_1d(spin_system,spec,parameters);

end

