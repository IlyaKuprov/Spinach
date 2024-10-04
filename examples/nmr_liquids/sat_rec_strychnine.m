% 1H saturation-recovery experiment on strychnine at 250 MHz. 
%
% Calculation time: minutes
%
% Zak El-Machachi
% i.kuprov@soton.ac.uk

function sat_rec_strychnine()

% Read spin system properties 
[sys,inter]=strychnine({'1H'});

% Magnetic induction
sys.magnet=5.9;

% Maximum distance to consider
sys.tols.prox_cutoff=5.0;

% Greedy parallelisation
sys.enable={'greedy'};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='kite';
inter.tau_c={200e-12};
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.offset=1250;
parameters.sweep=2500;
parameters.npoints=4096;
parameters.max_delay=0.5;
parameters.n_delays=10;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fids=liquid(spin_system,@sat_rec,parameters,'nmr');

% Apodisation
fids=apodisation(spin_system,fids,{{'exp',6},{}});

% Fourier transform
spectra=fftshift(fft(fids,[],1));

% Plotting
figure(); plot_1d(spin_system,real(spectra),parameters);

end

