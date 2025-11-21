% 13C{1H} HOESY spectrum of strychnine with natural content 
% of 13C isotope. 
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% zem1g16@soton.ac.uk

function hoesy_strychnine()

% Strychnine spin system
[sys,inter]=strychnine({'1H','13C'});

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.connectivity='scalar_couplings';
bas.space_level=3;
bas.level=4;

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
parameters.sweep=[6000 18000];
parameters.offset=[3000 12000];
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
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'positive');

end

