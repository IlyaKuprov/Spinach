% PANSY-COSY spectra of camphor with natural content of 13C isotope. 
% Coordinates, shieldings, and J-couplings computed with DFT.
%
% Calculation time: seconds
%
% Andrew Porter
% i.kuprov@soton.ac.uk

function pansy_dual_ch()

% Spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'),...
            {{'H','1H'},{'C','13C'}},[31.5 189.2],options); 
         
% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);

% Sequence parameters
parameters.sweep=[1800 9000];
parameters.offset=[900 4500];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'1H','13C'};
parameters.axis_units='ppm';

% Generate isotopomers
subsystems=dilute(spin_system,'13C');

% Preallocate the answer
spec_a=zeros(parameters.zerofill(1),parameters.zerofill(1),'like',1i);
spec_b=zeros(parameters.zerofill(2),parameters.zerofill(1),'like',1i);

% Loop over isotopomers
parfor n=1:numel(subsystems)
    
    % Build the basis
    subsystem=basis(subsystems{n},bas);
    
    % Simulation
    fid=liquid(subsystem,@pansy_cosy,parameters,'nmr');

    % Apodization
    fid.aa=apodization(fid.aa,'sqcosbell-2d');
    fid.ab=apodization(fid.ab,'sqcosbell-2d');

    % Fourier transforms
    spec_a=spec_a+fftshift(fft2(fid.aa,parameters.zerofill(1),...
                                       parameters.zerofill(1)));
    spec_b=spec_b+fftshift(fft2(fid.ab,parameters.zerofill(2),...
                                       parameters.zerofill(1)));

end

% Plotting - heteronuclear side
figure(); scale_figure([2.0 1.5]); subplot(1,2,2);
plot_2d(spin_system,abs(spec_b),parameters,20,...
        [0.05 0.25 0.05 0.25],2,256,6,'positive');

% Plotting - homonuclear side
subplot(1,2,1);
parameters.sweep=parameters.sweep(1);
parameters.offset=parameters.offset(1);
parameters.zerofill=parameters.zerofill(1);
parameters.spins=parameters.spins(1);
plot_2d(spin_system,abs(spec_a),parameters,20,...
        [0.05 0.25 0.05 0.25],2,256,6,'positive');

end

