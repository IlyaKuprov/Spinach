% 1H NMR spectrum of strychnine including an accurate model of
% line widths via Redfield superoperator.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function pa_strychnine()

% Read spin system properties 
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='kite';
inter.tau_c={200e-12};

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=2800;
parameters.sweep=6500;
parameters.npoints=8192;
parameters.zerofill=65536;
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

