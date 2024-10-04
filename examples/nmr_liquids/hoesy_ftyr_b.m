% (19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not 
% the right way to run this sequence in proteins because aro-
% matic 19F T2 is short, but 19F is being phase-encoded.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk

function hoesy_ftyr_b()

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
parameters.sweep=[2500 4000];
parameters.offset=[-70000 3000];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'19F','1H'};
parameters.decouple_f2={'19F'};
parameters.decouple_f1={'1H'};
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

