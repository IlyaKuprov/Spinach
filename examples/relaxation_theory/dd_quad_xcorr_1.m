% Complete Bloch-Redfield-Wangsness relaxation superoperator in a system 
% with a quadrupolar coupling and a dipole coupling. Spinach relaxation 
% theory module automatically accounts for all cross-correlations (dipole-
% quadrupole cross-correlation is present in this case). Dipolar couplings
% are computed from Cartesian coordinates of the two spins.
%
% Calculation time: seconds
%
% i.kuprov@soton.ac.uk

function dd_quad_xcorr_1()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','14N'};
inter.coupling.scalar={0 50; 50 0};
inter.coupling.eigs{2,2}=[1e4 1e4 -2e4];
inter.coupling.euler{2,2}=[0 0 0];
inter.coordinates={[0 0 0]
                   [0 0 1.02]};

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={1e-9};

% Basis specification
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=500;
parameters.npoints=128;
parameters.zerofill=512;
parameters.axis_units='Hz';

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'exp-1d',6);

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

