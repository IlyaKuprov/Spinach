% Methyl radical simulation, Gaussian import. The uncommon signal
% intensity pattern comes from g-HFC cross-correlation.
%
% luca.sementa@pi.ipcf.cnr.it

function gaussian_import_example() 

% System properties (vacuum DFT calculation)
options.no_xyz=1;
[sys,inter]=g2spinach(gparse('gaussian_methyl_radical.out'),...
                     {{'E','E'},{'H','1H'}},[0 0],options);

% Magnet induction
sys.magnet=0.339;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={5e-10};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'E'};
parameters.rho0=state(spin_system,'L+','E');
parameters.coil=state(spin_system,'L+','E');
parameters.decouple={};
parameters.sweep=5e8;
parameters.npoints=512;
parameters.zerofill=1024;
parameters.axis_units='mT';
parameters.derivative=1;
parameters.invert_axis=0;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'esr');

% Apodisation
fid=apodisation(spin_system,fid,{{'none'}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end
