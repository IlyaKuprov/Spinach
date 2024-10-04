% X-band pulse-acquire FFT ESR spectrum of methyl radical. Simple
% common line width is used as a relaxation model. Set to reprodu-
% ce Figure 4 from the paper by Zhitnikov and Dmitriev:
%
%        http://dx.doi.org/10.1051/0004-6361:20020268
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function pulse_acquire_methyl()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% Read the spin system (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/methyl.log'),...
                     {{'E','E'},{'H','1H'}},[0 0],options);
% Magnet induction
sys.magnet=0.33;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;
bas.longitudinals={'1H'};

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=2.5e7;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=5e8;
parameters.npoints=256;
parameters.zerofill=1024;
parameters.axis_units='GHz-labframe';
parameters.derivative=1;
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'none'}});

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

