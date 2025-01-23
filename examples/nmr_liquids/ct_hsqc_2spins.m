% CT HSQC spectrum of 2 spin system
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% mrw1g16@soton.ac.uk

function ct_hsqc_2spins()

% Spin system
sys.isotopes={'13C','1H'};

% Magnet field
sys.magnet=5.9;
inter.zeeman.scalar={50.00 3.00};
inter.coupling.scalar{1,2}=140.0; 
inter.coupling.scalar{2,2}=0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Sequence parameters
parameters.J=140;
parameters.sweep=[2500 950];
parameters.offset=[3000 600];
parameters.npoints=[128 128];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
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
    fid=liquid(subsystem,@ct_hsqc,parameters,'nmr');
    
    % Apodisation
    fid.pos=apodisation(spin_system,fid.pos,{{'cos'},{'cos'}});
    fid.neg=apodisation(spin_system,fid.neg,{{'cos'},{'cos'}});
    
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
plot_2d(spin_system,abs(spectrum),parameters,...
        20,[0.05 1.0 0.05 1.0],2,256,6,'positive');

end

