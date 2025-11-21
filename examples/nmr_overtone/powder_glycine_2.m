% Overtone detection 14N powder NMR spectrum of glycine, computed using
% Fokker-Planck formalism. Glycine quadrupolar tensor data comes from
% the paper by O'Dell and Ratcliffe:
%
%            http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% This simulation demonstrates that the spin state that gives rise to
% the overtone signal in a static sample is the T2,-2 coherence.
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function powder_glycine_2()

% System specification, TODO: add CSA
sys.magnet=14.1; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);
inter.zeeman.scalar{1}=32.4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=500;

% Algorithmic options
sys.disable={'krylov','trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Spectrum setup
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=[0e3 15e3];
parameters.npoints=256;
parameters.zerofill=256;
parameters.rho0=state(spin_system,'T2,-2','14N');
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*state(spin_system,'Lx','14N');
parameters.spins={'14N'};
parameters.axis_units='kHz';
parameters.verbose=0;

% Simulation
spectrum=powder(spin_system,@overtone_a,parameters,'qnmr');

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

