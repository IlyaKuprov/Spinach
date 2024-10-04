% 1H-13C RESPIRATION-CP experiment in the doubly rotating frame.
% Magic angle spinning simulation using Fokker-Planck formalism.
%
% Calculation time: seconds
%
% venkata-subbarao.redrouthu@uni-konstanz.de
% guinevere.mathies@uni-konstanz.de
% i.kuprov@soton.ac.uk

function cp_respiration()

% Magnet field
sys.magnet=11.7;

% System specification
sys.isotopes={'1H','13C'};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 2.00]};

% Formalism and basis
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.rate=20000;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.max_rank=8;
parameters.spins={'1H','13C'};
parameters.rho0=state(spin_system,'Lx','1H'); 
parameters.coil=state(spin_system,'L+','13C');
parameters.grid='rep_2ang_100pts_sph';
parameters.npoints=512;
parameters.zerofill=16384;
parameters.axis_units='kHz';
parameters.sweep=40000;
parameters.nloops=16;
parameters.theta=pi/20;

% Simulation
fid=singlerot(spin_system,@respiration,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',5}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); parameters.spins={'13C'};
plot_1d(spin_system,real(spectrum),parameters);

end

