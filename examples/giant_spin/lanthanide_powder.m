% Powder spectrum of Gd(III) with ZFS up to 3rd spherical rank
% using the giant spin Hamiltonian formalism.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function lanthanide_powder()

% Spin system properties
sys.isotopes={'E8'};
inter.zeeman.scalar={1.9918};
inter.giant.coeff={{[0 0 0],[0 0 -4.65e8 0 0],[1e7 0 0 2e7 0 0 1e7]}};
inter.giant.euler={{[0 0 0],[0 0 0],[0 0 0]}};

% Magnet field
sys.magnet=9.40;

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Disable trajectory-level SSR algorithms
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E8'};
parameters.rho0=state(spin_system,'L+','E8');
parameters.coil=state(spin_system,'L+','E8');
parameters.decouple={};
parameters.offset=1.5e9;
parameters.sweep=6e9;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='Gauss';
parameters.grid='rep_2ang_6400pts_sph';
parameters.rframes={{'E8',2}};
parameters.derivative=0;
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

