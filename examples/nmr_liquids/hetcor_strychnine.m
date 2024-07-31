% HETCOR spectrum of strychnine with natural content of 13C isotope.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function hetcor_strychnine()

% Read the spin system properties
[sys,inter]=strychnine({'1H','13C'});

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
parameters.J=140;
parameters.sweep=[3000 10000];
parameters.offset=[1000 4000];
parameters.npoints=[256 256];
parameters.zerofill=[512 512];
parameters.spins={'1H','13C'};
parameters.decouple={'1H'};
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
    fid=liquid(subsystem,@hetcor,parameters,'nmr');
    
    % Apodization
    fid=apodization(fid,'cosbell-2d');
    
    % Fourier transform
    spectrum=spectrum+fftshift(fft2(fid,parameters.zerofill(2),...
                                        parameters.zerofill(1)));
    
end

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.05 0.5 0.05 0.5],2,256,6,'positive');

end

