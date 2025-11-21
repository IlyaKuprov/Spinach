% X-band pulse-acquire FFT ESR spectrum of parafluorotoluene
% radical, simulated using explicit time-domain propagation
% including Redfield relaxation superoperator.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function relaxation_parafluorotoluene()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% Read the spin system (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/parafluorotoluene.log'),...
                      {{'E','E'},{'H','1H'},{'F','19F'}},[0 0 0],options);

% Ignore small HFC anisotropies
sys.tols.inter_cutoff=1e5;

% Magnet field
sys.magnet=0.33;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.zero_quantum={[1 2 3 5]};

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={1e-10};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=3e8;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='GHz-labframe';
parameters.derivative=1;
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'none'}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

