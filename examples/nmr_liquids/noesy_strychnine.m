% NOESY spectrum of strychnine.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de

function noesy_strychnine()

% Spin system properties
[sys,inter]=strychnine({'1H'});

% Magnet field
sys.magnet=5.9;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=3;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='IME';
inter.temperature=298;
inter.rlx_keep='kite';
inter.tau_c={200e-12};

% Algorithmic options
sys.enable={'greedy'};
sys.disable={'krylov'};
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.5;
parameters.offset=1200;
parameters.sweep=[2500 2500];
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.needs={'rho_eq'};

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

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
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.01 0.1 0.01 0.1],2,256,6,'both');

end

