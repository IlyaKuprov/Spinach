% Basic two-transmon system with Duffing model inter-
% actions and coupling: Coherence transfer from transmon 1 to transmon 2.
% System Hamiltonian and paramaters sourced from: https://doi.org/10.1088/2632-2153/ad4f4d
%
% c.musselwhite@soton.ac.uk
%

function two_transmons_RL()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','T3'};

% Rotating frame transmon parameters
inter.modes.frqs={-86.6e6 0};
inter.modes.anharms={-310.5e6 -313.9e6};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=2.2e6;

% Formalism and basis
bas.formalism='zeeman-liouv';
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

% Unit fidelity is Sorensen bound ('zeeman-hilb')
%rho_targ=rho_targ/sorensen(rho_init,rho_targ);

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={C_A,C_B};                      % Controls
control.off_ops={O_A,O_B};                        % Offset operator
control.offsets={linspace(-200,200,5)...
                 linspace(-200,200,5)};
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*[100e6 150e6 200e6 250e6 300e6];       % Pulse power ensemble
control.pulse_dt=1e-9*ones(1,300);                % Slice durations
control.penalties={'NS'};                         % Penalties
control.p_weights=0.1;                            % Penalty weights
control.method='goodwin';                         % Optimisation method
control.max_iter=5;                               % Termination condition
control.parallel='ensemble';                      % Parallelisation mode

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess - random
pulse=(1/8)*randn(2,300);

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

end

