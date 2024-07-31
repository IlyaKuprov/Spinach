% CLIP-HSQC spectrum of camphor with natural content of 13C isotope. 
% Coordinates, shielding anisotropies and J-couplings computed with
% DFT, isotropic chemical shifts taken from experimental data.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function clip_hsqc_camphor()

% Spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=0;
[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'),...
                 {{'H','1H'},{'C','13C'}},[31.8 182.1],options);          
% Magnet field
sys.magnet=14.1;
                     
% Isotropic shift components come from the experiment
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1:26,...
                              [ 29.70 26.80 44.20 59.70 42.80 ...
                               218.1  49.70 17.80 17.30  7.80 ...
                                 1.35  1.67  1.97  1.31  1.96 ...
                                 0.85  0.85  0.85  0.99  0.99 ...
                                 0.99  0.90  0.90  0.90  2.33 1.76]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Sequence parameters
parameters.J=140;
parameters.sweep=[8000 1500];
parameters.offset=[4000 1000];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
parameters.axis_units='ppm';

% Create the spin system structure
spin_system=create(sys,inter);

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
    fid=liquid(subsystem,@clip_hsqc,parameters,'nmr');
    
    % Apodization
    fid.pos=apodization(fid.pos,'sqcosbell-2d');
    fid.neg=apodization(fid.neg,'sqcosbell-2d');
    
    % F2 Fourier transform
    f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1);
    f1_neg=fftshift(fft(fid.neg,parameters.zerofill(2),1),1);
    
    % Form States signal
    fid=f1_pos+conj(f1_neg);
    
    % F1 Fourier transform
    spectrum=spectrum+fftshift(fft(fid,parameters.zerofill(1),2),2);
    
end

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 0.5 0.05 0.5],2,256,6,'negative');

end

