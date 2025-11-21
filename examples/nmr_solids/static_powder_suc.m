% 13C NMR spectrum of static sucrose powder, protons are assumed
% to be decoupled.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function static_powder_suc()

% Spin system properties (PCM DFT calculation)
[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'),...
                                        {{'C','13C'}},182.1,[]);
% Magnet field
sys.magnet=14.1;

% Disable trajectory-level SSR algorithms 
sys.disable={'trajlevel'};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.projections=+1;
bas.level=3;

% Algorithmic options
sys.tols.inter_cutoff=5.0;
sys.tols.prox_cutoff=4.0;
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment setup
parameters.sweep=5e4;
parameters.npoints=128;
parameters.zerofill=512;
parameters.offset=15000;
parameters.spins={'13C'};
parameters.decouple={};
parameters.axis_units='ppm';
parameters.invert_axis=1;
parameters.grid='rep_2ang_800pts_sph';
parameters.rho0=state(spin_system,'L+','13C','cheap');
parameters.coil=state(spin_system,'L+','13C','cheap');
parameters.verbose=0;

% Simulation
fid=powder(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

