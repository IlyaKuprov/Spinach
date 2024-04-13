% Optimal control pulse optimisation for state-to-state transfer 
% in a quadrupolar 14N spin at a fixed orientation and power level
%
% Two pulses are optimised cooperatively, so that the sum of the
% outcomes only contains the target state, and no impurities.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk

function state_transfer_coop()

% Magnet field
sys.magnet=14.1; 

% Isotopes
sys.isotopes={'14N'}; 

% Glycine NQI, random orientation
euler_angles=[1.0 2.0 3.0];
inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,euler_angles);

% Glycine 14N chemical shift
inter.zeeman.scalar{1}=32.4;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none'; 

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,'T1,0','14N');
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=state(spin_system,'T2,0','14N');
rho_targ=rho_targ/norm(full(rho_targ),2);

% spin_system assumptions
spin_system=assume(spin_system,'qnmr');

% Get the drift Hamiltonian
[Iso,Q]=hamiltonian(spin_system);
H=Iso+orientation(Q,[1 2 3]);
C=carrier(spin_system,'14N');
H=rotframe(spin_system,C,H,'14N',2);

% Get the control operators
Lx=operator(spin_system,'Lx','14N');
Ly=operator(spin_system,'Ly','14N');

% Define control parameters
control.drifts={{H}};               % Drift
control.operators={Lx,Ly};          % Controls
control.rho_init={rho_init};        % Starting state
control.rho_targ={rho_targ};        % Destination state
control.pwr_levels=2*pi*50e3;       % Pulse power
control.pulse_dt=10e-8*ones(1,100); % Slice durations
control.method='lbfgs';             % Optimisation method
control.amplitudes=ones(1,100);     % Amplitude profile
control.max_iter=100;               % Termination condition

% Plotting options
control.plotting={'coherence_order','phi_controls'};

% Initial guess
guess=randn(2,100);

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Optimisation
[~,data]=fminnewton(spin_system,@grape_coop,guess);

% Pull out final states
outcome_a=data.traj_data{1}{1}.forward(:,end);
outcome_b=data.traj_data{2}{1}.forward(:,end);

% Print diagnostics
disp('   Outcome A          Outcome B          (A+B)/2            Target');
disp(full([outcome_a outcome_b (outcome_a+outcome_b)/2 rho_targ]));

end

