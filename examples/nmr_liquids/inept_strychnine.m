% INEPT experiment on strychnine.
%
% Calculation time: minutes
%
% Andrew Porter, Ilya Kuprov

function inept_strychnine()

% Read the spin system properties 
[sys,inter]=strychnine({'1H','13C'});

% Magnet field
sys.magnet=5.9;

% Temperature
inter.temperature=298;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

%Sequence parameters
parameters.sweep=10000;
parameters.offset=[5000 0];
parameters.npoints=2048;
parameters.zerofill=8196;
parameters.spins={'13C','1H'};
parameters.J=150;
parameters.beta=3*pi/4;

% Create the spin system structure
spin_system=create(sys,inter);

% Generate isotopomers
subsystems=dilute(spin_system,'13C');

% Preallocate the answer
spectrum=zeros([parameters.zerofill 1],'like',1i);

% Loop over isotopomers
parfor n=1:numel(subsystems)
     
    % Build the basis
    subsystem=basis(subsystems{n},bas);
    
    % Simulation
    fid=liquid(subsystem,@inept,parameters,'nmr');
    
    % Apodization
    fid=apodization(fid,'exp-1d',6);
    
    % Fourier transform
    spectrum=spectrum+fftshift(fft(fid,parameters.zerofill));
    
end

%Plotting
figure(); parameters.spins={'13C'};
parameters.offset=parameters.offset(1);
parameters.axis_units='ppm';
plot_1d(spin_system,imag(spectrum),parameters);

end

