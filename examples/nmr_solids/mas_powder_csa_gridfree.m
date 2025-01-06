% Powder magic angle spinning spectrum of a pair of 
% anisotropically shielded proton spins using grid-
% free Fokker-Planck equation formalism.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function mas_powder_csa_gridfree()

% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5};
inter.zeeman.euler={[0 0 0],[0 0 0]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
bas.projections=+1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Parameters
parameters.axis=[1 1 1];
parameters.sweep=2e4;
parameters.npoints=512;
parameters.zerofill=4096;
parameters.offset=0;
parameters.spins={'1H'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.rate=500;
parameters.max_rank=17;
parameters.verbose=0;

% Simulation
fid=gridfree(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

