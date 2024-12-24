% Cross-polarization experiment between protons and 14N overtone
% transition in glycine under MAS. Glycine quadrupolar tensor da-
% ta comes from the paper by O'Dell and Ratcliffe:
%
%         http://dx.doi.org/10.1016/j.cplett.2011.08.030
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il
% m.carravetta@soton.ac.uk
% m.concistre@soton.ac.uk

function cpmas_glycine_simple()

% System specification
sys.magnet=14.1; sys.isotopes={'14N','1H'};
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0]);
inter.coupling.matrix{2,2}=[];
inter.zeeman.scalar={32.4  0};
inter.coordinates={[0    0    0]
                   [0    0    1.00]};
                  
% Relaxation theory
inter.relaxation={'damp'};
inter.rlx_keep='diagonal';
inter.equilibrium='zero';
inter.damp_rate=300;

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

% Spectrum setup
parameters.max_rank=7;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=-19840;
parameters.grid='rep_2ang_6400pts_sph';
parameters.sweep=[4.4e4 5.2e4];
parameters.npoints=256;
parameters.zerofill=256;
parameters.rho0=cos(theta)*state(spin_system,'Lz','1H')+...
                sin(theta)*(state(spin_system,'L+','1H')+...
                            state(spin_system,'L-','1H'))/2;
parameters.coil=cos(theta)*state(spin_system,'Lz','14N')+...
                sin(theta)*(state(spin_system,'L+','14N')+...
                            state(spin_system,'L-','14N'))/2;
parameters.Nx=cos(theta)*operator(spin_system,'Lz','14N')+...
              sin(theta)*(operator(spin_system,'L+','14N')+...
                          operator(spin_system,'L-','14N'))/2;
parameters.Hx=cos(theta)*operator(spin_system,'Lz','1H')+...
              sin(theta)*(operator(spin_system,'L+','1H')+...
                          operator(spin_system,'L-','1H'))/2;
parameters.spins={'14N'};
parameters.method='average';
parameters.axis_units='kHz';
parameters.rf_pwr=2*pi*[55.0e3 35.1e3]/sin(theta);
parameters.rf_dur=1e-4;
parameters.rf_frq=48e3;
parameters.verbose=0;

% Simulation
spectrum=singlerot(spin_system,@overtone_cp,parameters,'qnmr');

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

