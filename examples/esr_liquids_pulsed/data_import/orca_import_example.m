% Methyl radical simulation, ORCA import. The uncommon signal
% intensity pattern comes from g-HFC cross-correlation.
%
% luca.sementa@pi.ipcf.cnr.it

function orca_import_example() 

% System properties (vacuum DFT calculation)
props=oparse('orca_methyl_radical.out');

% Isotopes
sys.isotopes={'E','1H','1H','1H'};

% Zeeman interactions
inter.zeeman.matrix=cell(1,4);
inter.zeeman.matrix{1}=props.g_tensor.matrix;

% Hyperfine couplings
inter.coupling.matrix=cell(4,4);
inter.coupling.matrix{1,2}=1e6*gauss2mhz(props.hfc.full.matrix{2});
inter.coupling.matrix{1,3}=1e6*gauss2mhz(props.hfc.full.matrix{3});
inter.coupling.matrix{1,4}=1e6*gauss2mhz(props.hfc.full.matrix{4});

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
