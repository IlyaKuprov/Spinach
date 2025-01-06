% Powder magic angle spinning CN2D experiment (rotor-synchronized de-
% tection) on a 14N-1H spin pair using 1D Fokker-Planck equation and
% a spherical grid. The calculation accounts for the second-order qu-
% adrupolar shift and lineshape.
%
% Calculation time: hours on CPU, minutes with a Tesla V100 GPU.
%
% ilya.kuprov@weizmann.ac.il
% p.t.williamson@soton.ac.uk
% j.jarvis@soton.ac.uk

function hmqc_mas_dq()

% System specification
sys.magnet=19.96; 
sys.isotopes={'14N','1H'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);
inter.coupling.matrix{2,2}=[];
inter.zeeman.scalar={32.4  5};
inter.coordinates={[0.00 0.00 0.00]
                   [0.00 0.00 1.00]};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Use GPU if present
sys.enable={'gpu'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.rate=125000;
parameters.axis=[1 1 1];
parameters.max_rank=16;
parameters.grid='rep_2ang_200pts_sph';
parameters.sweep=[parameters.rate/4 20000];
parameters.npoints=[256 128];
parameters.zerofill=[1024 512];
parameters.spins={'14N','1H'};
parameters.rframes={{'14N',3}};
parameters.axis_units='ppm';
parameters.rho0=state(spin_system,'L+',parameters.spins{2},'cheap')+...
                state(spin_system,'L-',parameters.spins{2},'cheap');
parameters.coil=state(spin_system,'L+',parameters.spins{2},'cheap');
parameters.rf_pwr=40e3;
parameters.rf_dur=2e-3;

% Simulation
fid=singlerot(spin_system,@cn2d_dq,parameters,'qnmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1);
f1_sin=fftshift(fft(fid.sin,parameters.zerofill(2),1),1);

% Form States signal
f1_states=real(f1_cos)+1i*real(f1_sin);

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 0.5 0.05 0.5],2,256,6,'positive');

end

