% (1H) -> (19F) HOESY spectrum of fluorotyrosine, with the magneti-
% sation transfer direction picked so as to minimise the time that
% 19F spends in the transverse plane. This is the only way to run
% this sequence in proteins because aromatic 19F T2 is short.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il

function hoesy_ftyr_a()

% Read 3-fluorotyrosine DFT calculation
[sys,inter]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'),...
                            {{'H','1H'},{'F','19F'}},[31.82 192.97]);
         
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=3;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='IME';
inter.rlx_keep='kite';
inter.tau_c={10e-9};            % large protein
inter.temperature=298;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=5.0;
sys.tols.inter_cutoff=2.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.5;            % quite long
parameters.sweep=[4000 2500];
parameters.offset=[3000 -70000];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'1H','19F'};
parameters.decouple_f2={'1H'};
parameters.decouple_f1={'19F'};
parameters.axis_units='ppm';
parameters.needs={'rho_eq'};

% Simulation
fid=liquid(spin_system,@hoesy,parameters,'nmr');

% Apodisation
fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

% F2 Fourier transform
f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

% Form States signal
f1_states=f1_cos-1i*f1_sin;

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'negative');

end

