% Experimental NOESY spectrum of the Harvard RNA.
%
% Shunsuke Imai
% Scott Robson
% Gerhard Wagner
% Zenawi Welderufael
% Ilya Kuprov

function rna_noesy_expt()

% Magnet field
spin_system.inter.magnet=17.62;

% Sequence parameters
parameters.offset=3473;
parameters.sweep=[7500.00 7496.252];
parameters.zerofill=[1024 4096];
parameters.spins={'1H','1H'};
parameters.axis_units='ppm';

% Load data
load('rna_noesy_expt.mat','spec_expt')

% Plotting
spin_system.sys.disable={}; spin_system.sys.output=1; 
kfigure(); scale_figure([1.5 2.0]);
plot_2d(spin_system,-spec_expt,parameters,...
        20,[0.001 0.05 0.001 0.05],2,256,6,'positive');

end

