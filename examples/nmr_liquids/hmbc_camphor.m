% HMBC spectrum of camphor with natural content of 13C isotope. 
% Coordinates, shielding anisotropies and J-couplings computed 
% witt DFT.
%
% Calculation time: seconds
%
% Bud Macaulay
% i.kuprov@soton.ac.uk

function hmbc_camphor()

% Spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=0;
[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'),...
            {{'H','1H'},{'C','13C'}},[31.8-0.35 182.1+7.14],options); 
         
% Magnet field
sys.magnet=14.1;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Sequence parameters
parameters.J=140;
parameters.delta_b=60e-3;
parameters.sweep=[40000 1500];
parameters.offset=[18000 900];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
parameters.axis_units='Hz';

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
    fid=liquid(subsystem,@hmbc,parameters,'nmr');
    
    % Apodization
    fid=apodization(fid,'cosbell-2d');
    
    % Fourier transform
    spectrum=spectrum+fftshift(fft2(fid,parameters.zerofill(2),...
                                        parameters.zerofill(1)));
    
end

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.05 0.25 0.05 0.25],2,256,6,'positive');

end

