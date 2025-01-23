% Overtone 10B magic angle spinning NMR spectrum. The sample is
% spinning in the JEOL direction. Parameters from Nghia Duong
% and Yusuke Nishiyama, realistic RF power and pulse width.
%
% Calculation time: hours
%
% ilya.kuprov@weizmann.ac.il

function mas_boron_2()

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
inter.damp_rate=100;

% Algorithmic options
sys.disable={'trajlevel'};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Magic angle
theta=atan(sqrt(2));

% Sequence parameters
parameters.max_rank=12;
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.rate=70000;
parameters.grid='rep_2ang_200pts_sph';
parameters.sweep=[-141e3 -139e3];
parameters.npoints=256;
parameters.zerofill=256;
parameters.rho0=state(spin_system,'Lz','10B');
parameters.coil=cos(theta)*state(spin_system,'Lz','10B')+...
                sin(theta)*(state(spin_system,'L+','10B')+...
                            state(spin_system,'L-','10B'))/2;
parameters.Lx=cos(theta)*operator(spin_system,'Lz','10B')+...
              sin(theta)*(operator(spin_system,'L+','10B')+...
                          operator(spin_system,'L-','10B'))/2;
parameters.spins={'10B'};
parameters.axis_units='kHz';
parameters.rf_pwr=2*pi*50e3/sin(theta);
parameters.rf_dur=2e-3;
parameters.rf_frq=-140e3;
parameters.method='average';
parameters.verbose=1;

% Simulation
spectrum=singlerot(spin_system,@overtone_pa,parameters,'qnmr');

% Phasing
spectrum=exp(1i*1.45)*spectrum;

% Plotting
figure(); plot_1d(spin_system,real(spectrum),parameters);

end

