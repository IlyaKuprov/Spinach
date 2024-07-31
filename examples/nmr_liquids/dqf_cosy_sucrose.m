% DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT).
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function dqf_cosy_sucrose()

% Spin system properties (vacuum DFT calculation)
options.min_j=2.0; options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'),...
                                     {{'H','1H'}},31.8,options);
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

% Sequence parameters
parameters.offset=800;
parameters.sweep=1700;
parameters.npoints=[512 512];
parameters.zerofill=[2048 2048];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Simulation
fid=liquid(spin_system,@dqf_cosy,parameters,'nmr');

% Apodization
fid.cos=apodization(fid.cos,'cosbell-2d');
fid.sin=apodization(fid.sin,'cosbell-2d');

% F2 Fourier transform
f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1);
f1_sin=fftshift(fft(fid.sin,parameters.zerofill(2),1),1);

% Form States signal
f1_states=real(f1_cos)-1i*real(f1_sin);

% F1 Fourier transform
spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2);

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-real(spectrum),parameters,...
        20,[0.02 0.2 0.02 0.2],2,256,6,'both');

end

