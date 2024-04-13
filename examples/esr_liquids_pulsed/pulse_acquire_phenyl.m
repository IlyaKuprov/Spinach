% W-band pulse-acquire FFT ESR spectrum of phenyl radical. Simple
% fixed line width is used as a relaxation model.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function pulse_acquire_phenyl()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% Read the spin system properties (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/phenyl.log'),...
                     {{'E','E'},{'H','1H'}},[0 0],options);
% Magnet induction
sys.magnet=3.5;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.longitudinals={'1H'};
bas.projections=+1;

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=1e7;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E','cheap');
parameters.coil=state(spin_system,'L+','E','cheap');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=2e8;
parameters.npoints=512;
parameters.zerofill=1024;
parameters.axis_units='GHz-labframe';
parameters.derivative=1;
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodization
fid=apodization(fid,'none-1d');

% Perform Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plot the spectrum
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

