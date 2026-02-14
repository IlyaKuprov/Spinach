% Optimal control pulse optimisation for state-to-state transfer across a
% scalar coupling in a hydrofluorocarbon fragment spin system. The start-
% ing state is Z-magnetisation on 1H, the destination state is Z-magneti-
% sation on 19F. There are six control channels and a distribution over
% control powers.
%
% A keyhole condition is specified: the system must be in a two-spin or-
% der at point 20 in the pulse time grid. LBFGS quasi-Newton optimiser is
% used with a piecewise-linear control sequence.
%
% Calculation time: minutes.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il

function features_keyhole()

% Magnetic field
sys.magnet=9.4;

% Spin system
sys.isotopes={'1H','13C','19F'};

% Chemical shifts, ppm
inter.zeeman.scalar={0.0, 0.0, 0.0};

% Scalar couplings, Hz (literature values)
inter.coupling.scalar=cell(3);
inter.coupling.scalar{1,2}=140;
inter.coupling.scalar{2,3}=-160;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,{'Lz'},{1});
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=state(spin_system,{'Lz'},{3});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Control operators
LxH=operator(spin_system,'Lx','1H');
LyH=operator(spin_system,'Ly','1H');
LxC=operator(spin_system,'Lx','13C');
LyC=operator(spin_system,'Ly','13C');
LxF=operator(spin_system,'Lx','19F');
LyF=operator(spin_system,'Ly','19F');

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={LxH,LyH,LxC,LyC,LxF,LyF};      % Controls
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*linspace(0.8e3,1.2e3,5);  % Pulse powers, rad/s
control.pulse_dt=2e-4*ones(1,50);                 % Interval durations
control.penalties={'SNS'};                        % Penalty
control.p_weights=100;                            % Penalty weight
control.method='lbfgs';                           % Optimisation method
control.max_iter=200;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation

% Set the keyhole
control.keyholes=cell(1,50);
control.keyholes{20}=@(rho)correlation(spin_system,rho,2,'all');

% Plots during optimisation
control.plotting={'correlation_order','local_each_spin',...
                  'xy_controls','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=randn(6,50)/3;

% Run the optimal control procedure
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Denormalise and format the waveform
pulse=pulse*mean(control.pwr_levels);
pulse=mat2cell(pulse,[1 1 1 1 1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,H,control.operators,pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');
fidelity=real(rho_targ'*rho);
report(spin_system,['Re[<target|rho(T)>] = ' num2str(fidelity)]);

end

