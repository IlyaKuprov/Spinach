% Menthol NMR spectrum from Damien Jeannerat.
%
% Calculation time: minutes.
%
% damien.jeannerat@unige.ch

function pa_menthol()

% System and interaction specification
load('menthol.mat','sys','inter');

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.connectivity='scalar_couplings';
bas.approximation='IK-2';
bas.space_level=1;
bas.projections=+1;
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[4 5 6],[9 10 11],[12 13 14]};

% Algorithms
sys.enable={'greedy'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters - 1H
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');
parameters.decouple={};
parameters.offset=1000;
parameters.sweep=2000;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'gauss',10}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

