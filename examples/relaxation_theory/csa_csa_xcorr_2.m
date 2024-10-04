% CSA-CSA cross-correlation in the 103Rh subsystem and its effect 
% on the widths of the three lines of the proton triplet.
%
% Calculation time: seconds.
% 
% fares@kofo.mpg.de

function csa_csa_xcorr_2()

% Magnet field
sys.magnet=11.75;

% Set the spin system
sys.isotopes={'1H','103Rh','103Rh'};
inter.zeeman.matrix={[6.9 0 0; 0 6.9 0; 0 0 6.9]...
                     [7250 0 0; 0 8000 0; 0 0 7250]...
                     [7250 0 0; 0 8000 0; 0 0 7250]};

% J-couplings
inter.coupling.scalar={0   4   4;
                       4   0   100; 
                       4   100 0};

% Relaxation theory
inter.relaxation={'redfield'};
inter.rlx_keep='secular';
inter.equilibrium='zero';
inter.tau_c={10e-9};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters - 1H
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=6.9*500;
parameters.sweep=50;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',20}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

