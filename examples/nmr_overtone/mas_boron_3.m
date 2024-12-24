% Overtone 10B magic angle spinning NMR spectrum. The sample is
% spinning in the JEOL direction. Parameters from Nghia Duong
% and Yusuke Nishiyama. An unphysically strong pulse is used 
% to obtain a panoramic spectrum. 
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function mas_boron_3()

% System specification
sys.magnet=16.4; sys.isotopes={'10B'};
inter.coupling.matrix{1,1}=eeqq2nqi(0.7e6,0.0,3,[0 0 0]);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=1000;

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.max_rank=12;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=70000;
parameters.grid='rep_2ang_200pts_oct';
parameters.sweep=[-200e3 200e3];
parameters.npoints=4096;
parameters.zerofill=4096;
parameters.rho0=state(spin_system,'Lz','10B');
parameters.coil=state(spin_system,'Lz','10B');
parameters.spins={'10B'};
parameters.axis_units='kHz';
parameters.verbose=1;

% Simulation
spectrum=singlerot(spin_system,@overtone_a,parameters,'qnmr');

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

