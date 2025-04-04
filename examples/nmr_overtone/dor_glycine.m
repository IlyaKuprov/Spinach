% Panoramic double rotation overtone 14N spectrum of glycine, 
% simulated as described in our paper (Figure 1B):
%
%             http://dx.doi.org/10.1039/C5CP03266K
%
% A short pulse with instrumentally inaccessible power is gi-
% ven to make the excitation pattern uniform. Glycine quadru-
% polar tensor data comes from O'Dell and Ratcliffe:
%
%        http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function dor_glycine()

% System specification
sys.magnet=14.1; sys.isotopes={'14N'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);

% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=100;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Algorithmic options
sys.disable={'krylov','trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Experiment setup
parameters.rate_outer=1425;
parameters.rate_inner=6950;
parameters.rank_outer=5;
parameters.rank_inner=5;
parameters.axis_outer=[sin(theta) 0 cos(theta)];                   % 54.74 degrees
parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))]; % 30.56 degrees
parameters.grid='rep_2ang_200pts_sph';
parameters.sweep=[-2e4 3e4];
parameters.npoints=1024;
parameters.zerofill=1024;
parameters.spins={'14N'};
parameters.rho0=state(spin_system,'Lz','14N');
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*state(spin_system,'Lx','14N');
parameters.Lx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*operator(spin_system,'Lx','14N');
parameters.method='average';
parameters.rf_pwr=2*pi*3.0e6;
parameters.rf_dur=1.0e-6;
parameters.rf_frq=-10e3;
parameters.axis_units='kHz';
parameters.verbose=0;

% Run the simulation
spectrum=doublerot(spin_system,@overtone_pa,parameters,'qnmr');

% Phasing
spectrum=exp(1i*1.49)*spectrum;

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

