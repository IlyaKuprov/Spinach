% 14N overtone DANTE spectrum of glycine, computed using Fokker-Planck
% formalism. Glycine quadrupolar tensor data comes from the paper by
% O'Dell and Ratcliffe:
%
%            http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% Calculation time: hours
%
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function dante_glycine()

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
parameters.max_rank=7;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=-19840;
parameters.grid='rep_2ang_1600pts_sph';
parameters.sweep=[-60e3 80e3];
parameters.npoints=2048;
parameters.zerofill=2048;
parameters.rho0=state(spin_system,'Lz','14N');
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*(state(spin_system,'L+','14N')+...
                            state(spin_system,'L-','14N'))/2;
parameters.Lx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*(operator(spin_system,'L+','14N')+...
                          operator(spin_system,'L-','14N'))/2;
parameters.spins={'14N'};
parameters.axis_units='kHz';
parameters.pulse_amp=2*pi*55e3/sin(theta);
parameters.pulse_dur=10e-6;
parameters.rf_frq=48e3;
parameters.method='average';
parameters.n_periods=4;
parameters.pulse_num=2;
parameters.verbose=0;

% Simulation
spectrum=singlerot(spin_system,@overtone_dante,parameters,'qnmr');

% Phasing
spectrum=exp(-1i*2.12)*spectrum;

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

