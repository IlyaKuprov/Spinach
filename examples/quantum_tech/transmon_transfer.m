% Basic two-transmon system with Duffing model interacti-
% ons and a flip-flop coupling; coherence transfer from
% transmon 1 to transmon 2. GRAPE optimisation with a
% distribution control powers and transmon offsets with
% a penalty on excess power.
%
% Calculation time: minutes.
%
% c.musselwhite@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function transmon_transfer()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3','T5'};

% Rotating frame transmon parameters
inter.modes.frqs={100e6 -200e6};
inter.modes.anharms={-10e6 -20e6};
inter.modes.exchange=cell(2,2);
inter.modes.exchange{1,2}=50e6;

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
control.isotopes={'T3','T5'};                          % Isotopes
control.channels=[1;2];                                % Channel map
control.drifts={{H}};                                  % Drift
control.operators={C_A,C_B};                           % Controls
control.off_ops={O_A,O_B};                             % Offset operators
control.offsets={linspace(-10e6,10e6,5)...
                 linspace(-10e6,10e6,5)};              % Offset distributions
control.rho_init={rho_init};                           % Starting state
control.rho_targ={rho_targ};                           % Destination state
control.pwr_levels=2*pi*[40e6 45e6 50e6 55e6 60e6]*5;  % Pulse power ensemble
control.pulse_dt=0.25e-9*ones(1,200);                  % Slice durations
control.penalties={'SNSA'};                            % Penalties
control.p_weights=1.0;                                 % Penalty weights
control.method='rbfgs';                                % Optimisation method
control.max_iter=200;                                  % Termination condition

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess - random
pulse=(1/10)*randn(2,200);

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

end

