% Basic two-transmon system with Duffing model inter-
% actions and a flip-flop coupling.
%
% ilya.kuprov@weizmann.ac.il

function two_transmons()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','T5'};

% Rotating frame transmon parameters
inter.modes.frqs={100 -200};
inter.modes.anharms={-10 -20};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=50;

% Formalism and basis
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Get elementary operators
CrA=operator(spin_system,'C',1);
AnA=operator(spin_system,'A',1);
CrB=operator(spin_system,'C',2);
AnB=operator(spin_system,'A',2);

% Build control operators
C_A=(CrA+AnA)/2; C_B=(CrB+AnB)/2;

% Build offset operators
O_A=operator(spin_system,'N',1);
O_B=operator(spin_system,'N',2);

% Build source and destination states
rho_init=state(spin_system,{'C','BL1'},{1 2})+...
         state(spin_system,{'A','BL1'},{1 2});
rho_targ=state(spin_system,{'BL1','C'},{1 2})+...
         state(spin_system,{'BL1','A'},{1 2});
rho_init=rho_init/norm(rho_init,'fro');
rho_targ=rho_targ/norm(rho_targ,'fro');

% Clean up rounding errors
rho_init=(rho_init+rho_init')/2;
rho_targ=(rho_targ+rho_targ')/2;

% Unit fidelity is Sorensen bound
rho_targ=rho_targ/sorensen(rho_init,rho_targ);

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={C_A,C_B};                      % Controls
control.off_ops={O_A,O_B};                        % Offset operator
control.offsets={linspace(-10,10,5)...
                 linspace(-10,10,5)};             % Offset distribution
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*[40 45 50 55 60]*5;       % Pulse power ensemble
control.pulse_dt=1e-3*ones(1,100);                % Slice durations
control.penalties={'NS'};                         % Penalties
control.p_weights=0.1;                            % Penalty weights
control.method='lbfgs';                           % Optimisation method
control.max_iter=100;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation mode

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess - random
pulse=(1/10)*randn(2,100);

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

end

