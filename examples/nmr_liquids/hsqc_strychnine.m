% HSQC spectrum of strychnine with natural content of 13C isotope.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function hsqc_strychnine()

% Spin system properties 
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
parameters.sweep=[10000 3000];
parameters.offset=[4000 1000];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
parameters.decouple_f1={'1H'};
parameters.decouple_f2={'13C'};
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
    fid=liquid(subsystem,@hsqc,parameters,'nmr');
    
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
        20,[0.05 1.0 0.05 1.0],2,256,6,'positive');

end

