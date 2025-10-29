% 13C 2D PDSD spectrum of a simple test spin system.
%
% Calculation time: minutes, much faster on GPU.
%
% ilya.kuprov@weizmann.ac.il
% guinevere.mathies@uni-konstanz.de

function pdsd_simple()

% Magnet field
sys.magnet=21.1356;

% Isotopes
sys.isotopes={'13C','1H','1H','13C'};

% Interactions (HCCH fragment)
inter.zeeman.scalar={20.0 5.0 2.0 35.0};
inter.coordinates={[-2.26  0.15  0.00],...
                   [-1.90  0.66 -0.87],...
                   [-0.67  0.88  1.26],...
                   [-1.74  0.88  1.26]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.enable={'prop_cache'}; % 'gpu'

% Create the spin system structure
spin_system=create(sys,inter);

% Build the basis
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.spins={'1H','13C'};
parameters.rate=10000;
parameters.tmix=10e-3;
parameters.max_rank=11;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='rep_2ang_100pts_sph';
parameters.npoints=[256 256];
parameters.zerofill=[1024 1024];
parameters.offset=[3150 6750];
parameters.sweep=10000;

% Simulation
fid=singlerot(spin_system,@pdsd,parameters,'nmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform and real part
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);
spectrum=real(spectrum);

% Plotting
figure(); scale_figure([1.5 2.0]);
parameters.spins={'13C'};
parameters.offset=parameters.offset(2);
plot_2d(spin_system,spectrum,parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'positive'); 

end

