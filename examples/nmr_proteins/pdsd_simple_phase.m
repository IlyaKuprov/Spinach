% 13C 2D PDSD spectrum of a simple text spin system.
%
% Calculation time: hours.
%
% ilya.kuprov@weizmann.ac.il

function pdsd_simple_phase()

% Magnet field
sys.magnet=21.1356;

% Isotopes
sys.isotopes={'13C','1H','1H','13C'};

% Interactions
inter.zeeman.scalar={20.0 5.0 2.0 35.0};
inter.coordinates={[0.00 0.00 -3.00],...
                   [0.00 0.00 -1.00],...
                   [0.00 0.00  1.00],...
                   [0.00 0.00  3.00]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.enable={'greedy','gpu'};

% Create the spin system structure
spin_system=create(sys,inter);

% Build the basis
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'1H','13C'};
parameters.rate=100000;
parameters.tmix=1e-3;
parameters.max_rank=3;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='rep_2ang_200pts_oct';
parameters.npoints=[64 64];
parameters.zerofill=[1024 1024];
parameters.offset=[3150 6750];
parameters.sweep=10000;
parameters.verbose=1;
parameters.serial=0;
parameters.axis_units='ppm';

% Simulation
fid=singlerot(spin_system,@pdsd_phase,parameters,'nmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
parameters.spins={'13C'};
parameters.offset=parameters.offset(2);
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'both'); 

end

