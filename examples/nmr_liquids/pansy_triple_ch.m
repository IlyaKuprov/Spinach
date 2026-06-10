% Triple-channel PANSY experiment on glycine with natural 
% content of 13C isotope. Coordinates, shieldings, and J-
% couplings computed with DFT.
%
% Calculation time: seconds
%
% Andrew Porter
% ilya.kuprov@weizmann.ac.il

function pansy_triple_ch()

% Read the spin system properties (vacuum DFT calculation)        
options.min_j=3.0; options.no_xyz=1;
[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'),...
            {{'H','1H'},{'C','13C'},{'N','15N'}},[31.5 189.2 400.3],options); 

% Magnet field
sys.magnet=5.9;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.prox_level=1;

% Spinach housekeeping
spin_system=create(sys,inter);

% Sequence parameters
parameters.sweep=[1700 4500 8000];
parameters.offset=[900 3200 4000];
parameters.npoints=[128 128 128];
parameters.zerofill=[512 512 512];
parameters.spins={'1H','13C','15N'};
parameters.axis_units='ppm';

% Preallocate the answer
spec_a=zeros(parameters.zerofill(1),parameters.zerofill(1),'like',1i);
spec_b=zeros(parameters.zerofill(2),parameters.zerofill(1),'like',1i);
spec_c=zeros(parameters.zerofill(3),parameters.zerofill(1),'like',1i);

% Generate isotopomers
subsys=dilute(spin_system,'13C');

% Loop over isotopomers
for n=1:numel(subsys)
    
    % Generate isotopomers
    subsubsys=dilute(subsys{n,:},'15N');
    
    % Loop over isotopomers
    for k=1:numel(subsubsys)
    
           % Build the basis
           subsystem=basis(subsubsys{k},bas);
    
           % Simulation
           fid=liquid(subsystem,@pansy_triple,parameters,'nmr');

           % Apodisation
           fid.aa=apodisation(spin_system,fid.aa,{{'sqcos'},{'sqcos'}});
           fid.ab=apodisation(spin_system,fid.ab,{{'sqcos'},{'sqcos'}});
           fid.ac=apodisation(spin_system,fid.ac,{{'sqcos'},{'sqcos'}});

           % Fourier transforms
           spec_a=spec_a+fftshift(fft2(fid.aa,parameters.zerofill(1),...
                                       parameters.zerofill(1)));
           spec_b=spec_b+fftshift(fft2(fid.ab,parameters.zerofill(2),...
                                       parameters.zerofill(1)));
           spec_c=spec_c+fftshift(fft2(fid.ac,parameters.zerofill(3),...
                                       parameters.zerofill(1)));

     end

end

% Store the original spectral axes
plot_offset=parameters.offset;
plot_sweep=parameters.sweep;
plot_zerofill=parameters.zerofill;
plot_spins=parameters.spins;

% Plotting - 1H-15N block
kfigure(); scale_figure([3.0 1.5]); subplot(1,3,3);
parameters.offset=[plot_offset(1) plot_offset(3)];
parameters.sweep=[plot_sweep(1) plot_sweep(3)];
parameters.zerofill=[plot_zerofill(1) plot_zerofill(3)];
parameters.spins=plot_spins([1 3]);
plot_2d(spin_system,real(spec_c),parameters,20,...
        [0.05 0.25 0.05 0.25],2,256,6,'both');

% Plotting - 1H-13C block
subplot(1,3,2);
parameters.offset=[plot_offset(1) plot_offset(2)];
parameters.sweep=[plot_sweep(1) plot_sweep(2)];
parameters.zerofill=[plot_zerofill(1) plot_zerofill(2)];
parameters.spins=plot_spins([1 2]);
plot_2d(spin_system,real(spec_b),parameters,20,...
        [0.05 0.25 0.05 0.25],2,256,6,'both');

% Plotting - 1H-1H block
subplot(1,3,1);
parameters.sweep=plot_sweep(1);
parameters.offset=plot_offset(1);
parameters.zerofill=plot_zerofill(1);
parameters.spins=plot_spins(1);
plot_2d(spin_system,real(spec_a),parameters,20,...
        [0.05 0.25 0.05 0.25],2,256,6,'both');

end

