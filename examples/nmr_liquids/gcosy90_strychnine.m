% Gradient-selected COSY spectrum of strychnine.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% gareth.charnock@oerc.ox.ac.uk
% ledwards@cbs.mpg.de

function gcosy90_strychnine()

% Read the spin system properties
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=5.9;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.angle=pi/2;
parameters.offset=1200;
parameters.sweep=2200;
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.g_amp=3;
parameters.g_dur=2e-3;
parameters.g_stab_del=2e-4;
parameters.s_len=1.5;
parameters.pathway='P+N';

% Simulation
fid=liquid(spin_system,@gcosy,parameters,'nmr');

% Apodisation
fid.pos=apodisation(spin_system,fid.pos,{{'cos'},{'cos'}});
fid.neg=apodisation(spin_system,fid.neg,{{'cos'},{'cos'}});

% F2 Fourier transform
f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);

% Form echo/anti-echo signal
fid=f1_pos+conj(f1_neg);

% F1 Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2);

% Plotting
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.02 0.2 0.02 0.2],2,256,6,'both');

end
