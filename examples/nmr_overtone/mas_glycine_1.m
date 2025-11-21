% Overtone detection 14N magic angle spinning NMR spectrum of glycine,
% computed using Fokker-Planck formalism. Glycine quadrupolar tensor
% data comes from the paper by O'Dell and Ratcliffe:
%
%            http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% Simulation parameters are set to reproduce Figure 3b from our paper
% on the subject:
%
%                http://dx.doi.org/10.1039/C4CP03994G
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function mas_glycine_1()

% System specification
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
inter.damp_rate=300;

% Algorithmic options
sys.disable={'krylov','trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Spectrum setup
parameters.max_rank=6;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=-19840;
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=[44e3 52e3];
parameters.npoints=256;
parameters.zerofill=256;
parameters.rho0=state(spin_system,'Lz','14N');
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*state(spin_system,'Lx','14N');
parameters.Lx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*operator(spin_system,'Lx','14N');
parameters.spins={'14N'};
parameters.axis_units='kHz';
parameters.rf_pwr=2*pi*55e3/sin(theta);
parameters.rf_dur=260e-6;
parameters.rf_frq=48e3;
parameters.method='average';
parameters.verbose=0;

% Simulation
spectrum=singlerot(spin_system,@overtone_pa,parameters,'qnmr');

% Phasing
spectrum=exp(1i*1.35)*spectrum;

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

