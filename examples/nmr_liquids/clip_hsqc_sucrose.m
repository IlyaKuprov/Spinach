% CLIP-HSQC spectrum of sucrose with natural content of 13C isotope. 
% Coordinates, shielding anisotropies and J-couplings computed with
% DFT, isotropic chemical shifts taken from experimental data.
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de

function clip_hsqc_sucrose()

% Read the spin system properties (vacuum DFT calculation)
options.min_j=3.0; options.no_xyz=0;
[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'),...
                 {{'H','1H'},{'C','13C'}},[31.8 182.1],options);

% Set the isotropic parts of shielding tensors to experimental values
spin_numbers=[1:19 24:30];
new_shifts=[ 94.5  73.4  74.9  71.5  74.7  62.4  63.6 ...
            106.0  78.7  76.3  83.7  64.7  5.49  3.63 ...
             3.83  3.54  3.90  3.90  3.90  3.75  3.75 ...
             4.29  4.12  3.96  3.90  3.90];
inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,spin_numbers,new_shifts);

% Set the field strength
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;

% Algorithmic options
sys.enable={'greedy'};
sys.tols.prox_cutoff=4.0;

% Sequence parameters
parameters.sweep=[8000 2000];
parameters.offset=[12000 2700];
parameters.npoints=[256 256];
parameters.zerofill=[512 512];
parameters.spins={'13C','1H'};
parameters.J=140;
parameters.axis_units='ppm';

% Create the spin system structure
spin_system=create(sys,inter);

% Remove fast exchanging and uncoupled spins from the simulation
spin_system=kill_spin(spin_system,[20,21,22,23,31,32,33,34]);

% Generate isotopomers
subsystems=dilute(spin_system,'13C');

% Preallocate the result
spectrum=zeros(parameters.zerofill,'like',1i);

% Run the CLIP-HSQC simulation
parfor n=1:numel(subsystems)
    
    % Build the basis
    subsystem=basis(subsystems{n},bas);
    
    % Run simulation
    fid=liquid(subsystem,@clip_hsqc,parameters,'nmr');
    
    % Apodization
    P_term=apodization(fid.pos,'sqcosbell-2d');
    N_term=apodization(fid.neg,'sqcosbell-2d');
    
    % F2 Fourier transform (directly detected dimension)
    f1_P=fftshift(fft(P_term,parameters.zerofill(2),1),1);
    f1_N=fftshift(fft(N_term,parameters.zerofill(2),1),1);
    
    % Form States signal
    f1_tot=f1_P+conj(f1_N);
    
    % F1 Fourier transform (indirectly detected dimension)
    spectrum=spectrum+fftshift(fft(f1_tot,parameters.zerofill(1),2),2);
    
end

% Plotting
figure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,real(spectrum),parameters,...
        20,[0.02 0.2 0.02 0.2],2,256,6,'negative');

end

