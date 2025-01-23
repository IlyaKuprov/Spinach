% Powder averaged W-band pulsed ESR spectrum of Gd(III) DOTA 
% complex. Ideal pulse with a large numerical powder grid is
% used, along with the numerical second-order rotating frame
% transformation.
%
% Calculation time: minutes
%
% corzilius@solidstatednp.com
% ilya.kuprov@weizmann.ac.il

function hpa_gd_dota_powder()

% Spin system properties
sys.isotopes={'E8'};
inter.zeeman.scalar={1.9918};
inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3;
inter.coupling.euler{1,1}=[0 0 0];

% Magnet field
sys.magnet=9.40;

% Basis set
bas.formalism='zeeman-liouv';
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
parameters.npoints=4096;
parameters.zerofill=16384;
parameters.axis_units='GHz-labframe';
parameters.grid='rep_2ang_12800pts_sph';
parameters.rframes={{'E8',2}};
parameters.derivative=0;
parameters.invert_axis=1;

% Simulation
fid=powder(spin_system,@acquire,parameters,'labframe');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

