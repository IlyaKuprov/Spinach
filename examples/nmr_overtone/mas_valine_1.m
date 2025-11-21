% Overtone detection 14N magic angle spinning NMR spectrum of 
% N-acetylvaline, computed using Fokker-Planck formalism. Va-
% line quadrupolar tensor data comes from our paper:
%
%            http://dx.doi.org/10.1039/c4cp03994g
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function mas_valine_1()

% System specification
sys.magnet=14.102; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(3.21e6,0.27,1,[0 0 0]);
inter.zeeman.eigs={[57.5 81.0 227.0]};
inter.zeeman.euler={[-90 -90 -17]*(pi/180)};

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=2000;

% Algorithmic options
sys.disable={'krylov','trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Spectrum setup
parameters.max_rank=8;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=-19840;
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=[75e3 100e3];
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
parameters.rf_dur=70e-6;
parameters.rf_frq=86e3;
parameters.method='average';
parameters.verbose=0;

% Simulation
spectrum=singlerot(spin_system,@overtone_pa,parameters,'qnmr');

% Phasing
spectrum=exp(1i*1.75)*spectrum;

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);

end

