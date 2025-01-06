% 1H NMR spectrum of styrene, the calculation is performed in Hilbert
% space to demonstrate parallel propagation described in:
%
%                 https://doi.org/10.1063/1.3679656
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il 

function pa_styrene()

% Isotopes
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'};

% Magnetic induction
sys.magnet=7.046;

% Chemical shifts
inter.zeeman.scalar={6.720  6.720  6.587  6.587  6.527  6.063   5.085  4.555};

% Scalar couplings
inter.coupling.scalar={ 0.0       1.9170    7.7884    0.6061    1.2447   -0.5347    0.0390    0.1582
                        0.0       0.0       0.6061    7.7884    1.2447   -0.5347    0.0390    0.1582
                        0.0       0.0       0.0       1.4184    7.4390    0.3693   -0.0416    0.0169
                        0.0       0.0       0.0       0.0       7.4390    0.3693   -0.0416    0.0169
                        0.0       0.0       0.0       0.0       0.0      -0.2279    0.2349    0.2905
                        0.0       0.0       0.0       0.0       0.0       0.0      17.6002   10.8954
                        0.0       0.0       0.0       0.0       0.0       0.0       0.0       1.0382
                        0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0};

% Basis set
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=1700;
parameters.sweep=1000;
parameters.npoints=16384;
parameters.zerofill=65536;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'gauss',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

