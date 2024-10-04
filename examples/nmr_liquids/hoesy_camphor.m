% 13C{1H} HOESY spectrum of camphor with natural content of 13C isotope. 
% Coordinates, shielding anisotropies and J-couplings computed with DFT.
%
% Calculation time: minutes
%
% Zak El-Machachi
% i.kuprov@soton.ac.uk

function hoesy_camphor()

% Spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=0;
[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'),...
            {{'H','1H'},{'C','13C'}},[31.5 189.2],options); 
         
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
inter.tau_c={50e-12};
inter.temperature=298;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=5.0;
sys.tols.inter_cutoff=2.0;

% Spinach housekeeping
spin_system=create(sys,inter);

% Sequence parameters
parameters.tmix=0.5;
parameters.sweep=[1800 9000];
parameters.offset=[900 4500];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'1H','13C'};
parameters.decouple_f2={'1H'};
parameters.decouple_f1={'13C'};
parameters.axis_units='ppm';
parameters.needs={'rho_eq'};

% Generate isotopomers
subsystems=dilute(spin_system,'13C');

% Preallocate the answer
spectrum=zeros(parameters.zerofill(2),...
               parameters.zerofill(1),'like',1i);

% Loop over isotopomers
parfor n=1:numel(subsystems)
    
    % Build the basis
    subsystem=basis(subsystems{n},bas);
    
    % Simulation
    fid=liquid(subsystem,@hoesy,parameters,'nmr');

    % Apodisation
    fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}});
    fid.sin=apodisation(spin_system,fid.sin,{{'sqcos'},{'sqcos'}});

    % F2 Fourier transform
    f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1));
    f1_sin=real(fftshift(fft(fid.sin,parameters.zerofill(2),1),1));

    % Form States signal
    f1_states=f1_cos-1i*f1_sin;

    % F1 Fourier transform
    spectrum=spectrum+fftshift(fft(f1_states,parameters.zerofill(1),2),2);

end

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'positive');

end

