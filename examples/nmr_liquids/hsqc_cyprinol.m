% HSQC spectrum of cyprinol with natural abundance of 13C isotope.
%
% Calculation time: seconds
%
% Bud Macaulay
% ilya.kuprov@weizmann.ac.il

function hsqc_cyprinol()

% Spin system
[sys,inter,bas]=cyprinol();

% Magnet field
sys.magnet=11.7;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;
sys.tols.inter_cutoff=5.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-1';
bas.level=3;
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Sequence parameters
parameters.J=150;
parameters.sweep=[12000 2500];
parameters.offset=[5000 1250];
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
    
    % Apodisation
    fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}});
    fid.neg=apodisation(spin_system,fid.neg,{{'sqcos'},{'sqcos'}});
    
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

