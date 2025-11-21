% Experimental HNCO spectrum of human ubiquitin.
%
% Donghan Lee (Max Planck Institute)
% Ilya Kuprov (University of Southampton)

function noesy_ubiquitin_expt()

% Magnet field
spin_system.inter.magnet=21.1356;

% Sequence parameters
parameters.offset=4250;
parameters.sweep=10815;
parameters.zerofill=[1024 1024];
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Load data
load('noesy_ubiquitin_expt.mat','spectrum')

% Plotting
spin_system.sys.disable={'colorbar'}; 
spin_system.sys.output=1; 
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,spectrum,parameters,...
        20,[1e-5 50e-5 1e-5 50e-5],2,256,6,'positive');

end

