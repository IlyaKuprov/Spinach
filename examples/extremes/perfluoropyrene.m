% X-band pulsed ESR spectrum of perfluoropyrene cation radical, computed 
% using brute force operator algebra in the full 4,194,304 - dimensional
% Liouville space. This is deliberate - a much faster calculation is, of
% course, possible with a restricted basis set.
% 
% This calculation requires at least 64GB of RAM and illustrates the per-
% formance of trajectory-level state space restriction in Spinach.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function perfluoropyrene()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% Read the spin system properties (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/perfluoropyrene_cation.log'),...
                                        {{'E','E'},{'F','19F'}},[0 0],options);
% Magnet induction
sys.magnet=0.33;

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=2e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

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
parameters.npoints=2048;
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

