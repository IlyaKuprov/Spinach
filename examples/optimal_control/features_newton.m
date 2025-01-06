% Optimal control pulse optimisation for state-to-state transfer across
% two scalar couplings in a hydrofluorocarbon fragment spin system. The
% starting state is Lz on 1H, the destination state is Lz on 19F. There
% are six control channels, the pulse is designed to be stable with res-
% pect to proton transmitter offset and pulse nutation frequency drop.
%
% The optimisation uses Newton-Raphson GRAPE algorithm described in:
%
%                  http://dx.doi.org/10.1063/1.4949534
%
% with point-by-point variation and a penalty on the waveform exceeding
% a user-specified power threshold. The initial guess is a random pulse.
%
% Calculation time: minutes.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il

function features_newton()

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

% Offset operator
LzH=operator(spin_system,'Lz','1H');

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={LxH,LyH,LxC,LyC,LxF,LyF};    % Controls
control.off_ops={LzH};                          % Offset operator 
control.offsets={linspace(-1e3,+1e3,5)};        % Offset range, Hz
control.rho_init={rho_init};                    % Starting state
control.rho_targ={rho_targ};                    % Destination state
control.pwr_levels=2*pi*[0.8 0.9 1.0]*1e3;      % Pulse powers, rad/s
control.pulse_dt=1e-4*ones(1,100);              % Slice durations
control.penalties={'NS'};                       % Penalty type
control.p_weights=0.01;                         % Penalty weight
control.method='newton';                        % Optimisation method
control.max_iter=50;                            % Termination condition
control.parallel='ensemble';                    % Parallelisation

% Plots during optimisation
control.plotting={'correlation_order','local_each_spin',...
                  'xy_controls','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=randn(6,100)/10;   

% Run the optimal control procedure
pulse=fminnewton(spin_system,@grape_xy,guess);

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

