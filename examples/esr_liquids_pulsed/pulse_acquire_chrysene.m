% W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in
% a non-viscous liquid. Simple common line width is used as a relaxation
% model. Symmetry treatment is performed using the full S2xS2xS2xS2xS2xS2
% group direct product.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk
% matthew.krzystyniak@oerc.ox.ac.uk

function pulse_acquire_chrysene()

% Ignore coordinate information (HFCs provided)
options.no_xyz=1;

% Read the spin system (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/chrysene_cation.log'),...
                            {{'E','E'},{'H','1H'}},[0 0],options);
% Magnet induction
sys.magnet=3.5;

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=1e6;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.longitudinals={'1H'};
bas.projections=+1;

% Symmetry
bas.sym_spins={[1 7],[2 8],[3 9],[4 10],[5 11],[6 12]};
bas.sym_group={'S2','S2','S2','S2','S2','S2'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set the sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.offset=-2e7;
parameters.sweep=1e8;
parameters.npoints=1024;
parameters.zerofill=4096;
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

