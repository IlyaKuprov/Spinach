% HSQC spectrum of sucrose with natural content of 13C isotope
% (magnetic parameters computed with DFT).
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il
% ledwards@cbs.mpg.de

function hsqc_sucrose()

% Spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'),...
                 {{'H','1H'},{'C','13C'}},[31.8 182.1],options);       

% Set the isotropic parts of shielding tensors to experimental values
spin_numbers=[1:19 24:30];
new_shifts=[ 94.5  73.4  74.9  71.5  74.7  62.4  63.6 ...
            106.0  78.7  76.3  83.7  64.7  5.49  3.63 ...
             3.83  3.54  3.90  3.90  3.90  3.75  3.75 ...
             4.29  4.12  3.96  3.90  3.90];
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,spin_numbers,new_shifts);

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
parameters.sweep=[3350 950];
parameters.offset=[5000 1100];
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
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.05 1.0 0.05 1.0],2,256,6,'positive');

end

