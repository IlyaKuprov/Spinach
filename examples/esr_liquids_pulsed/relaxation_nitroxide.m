% W-band pulse-acquire FFT ESR spectrum of a nitroxide radical, 
% using explicit time domain simulation with Redfield relaxati-
% on supeoperator.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function relaxation_nitroxide()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% Spin system properties (imported from a DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'),...
                           {{'E','E'},{'N','14N'}},[0 0],options);
% Magnet induction
sys.magnet=3.5;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% RElaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={5e-11};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=-2e8;
parameters.sweep=2e8;
parameters.npoints=512;
parameters.zerofill=1024;
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

