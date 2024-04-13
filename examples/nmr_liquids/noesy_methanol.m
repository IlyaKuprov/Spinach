% NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker,
% CSA tensors from DFT. Note the presence of cross-peaks between 13C
% doublet components.
%
% Calculation time: seconds
%
% tim.claridge@chem.ox.ac.uk
% i.kuprov@soton.ac.uk

function noesy_methanol()

% Spin system properties (vacuum DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/methanol.log'),...
                     {{'H','1H'},{'C','13C'}},[31.8 182.4],[]);

% Remove the OH proton
sys.isotopes(end)=[];
inter.coordinates(end)=[];
inter.zeeman.matrix(end)=[];

% Put all chemical shifts on resonance
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,[1 2 3 4],[0 0 0 0]);

% Magnet field
sys.magnet=14.1;

% Assign J-couplings
inter.coupling.scalar=cell(4,4);
inter.coupling.scalar{1,2}=141;
inter.coupling.scalar{1,3}=141;
inter.coupling.scalar{1,4}=141;
inter.coupling.scalar{2,3}=-11;
inter.coupling.scalar{3,4}=-11;
inter.coupling.scalar{2,4}=-11;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='IME';
inter.temperature=298;
inter.rlx_keep='kite';
inter.tau_c={50e-12};

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.tmix=0.5;
parameters.offset=0;
parameters.sweep=[300 300];
parameters.npoints=[256 256];
parameters.zerofill=[1024 1024];
parameters.spins={'1H'};
parameters.axis_units='ppm';
parameters.needs={'rho_eq'};

% Simulation
fid=liquid(spin_system,@noesy,parameters,'nmr');

% Apodization
fid.cos=apodization(fid.cos,'sqcosbell-2d');
fid.sin=apodization(fid.sin,'sqcosbell-2d');

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
        20,[0.05 0.5 0.05 0.5],2,256,6,'both');

end

