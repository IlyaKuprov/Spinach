% 1H NMR spectrum of sucrose (magnetic parameters read in from a DFT 
% calculation), including Redfield relaxation superoperator.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function pa_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=1.0;
[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'),...
                                     {{'H','1H'}},31.8,options);
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=2;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={1e-9};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=1800;
parameters.sweep=5000;
parameters.npoints=8192;
parameters.zerofill=65536;
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

