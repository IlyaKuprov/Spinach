% Basic single transmon system with Duffing model inter-
% actions, parameters and model from:
%
%         https://doi.org/10.1038/ncomms10628
%
% Frequencies are in scaled units of MHz and times in
% scaled units of microseconds. Ensemble GRAPE opti-
% misation with a distribution of control powers.
%
% Calculation time: minutes
%
% c.musselwhite@soton.ac.uk

function transmon_stirap()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3'};

% Rotating frame ladder detunings
inter.modes.frqs={10};
inter.modes.anharms={-20};

% Formalism and basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Pulse power ensemble
pwr_levels=2*pi*[30 35 40 45 50];

% Gaussian STIRAP pulse pair, normalised to the top power level
t=1e-3*linspace(-150,150,300);
ts=-90e-3; sigma=45e-3;
omega01=2*pi*43.4/max(pwr_levels);
omega12=2*pi*38.2/max(pwr_levels);
omega01_t=omega01*exp(-((t.^2)/(2*sigma^2)));
omega12_t=omega12*exp(-(((t-(ts/2)).^2)/(2*sigma^2)));

% Build control operators
H01=sparse([0 1 0; 1 0 0; 0 0 0]/2);
H12=sparse([0 0 0; 0 0 1; 0 1 0]/2);

% Build source and destination states
rho_init=state(spin_system,'BL1',1);
rho_targ=state(spin_system,'BL3',1);
rho_init=rho_init/norm(rho_init,'fro');
rho_targ=rho_targ/norm(rho_targ,'fro');

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={hilb2liouv(H01,'comm'),...
                   hilb2liouv(H12,'comm')};       % Controls
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=pwr_levels;                    % Pulse power ensemble
control.pulse_dt=(t(2)-t(1))*ones(1,300);         % Slice durations
control.penalties={'NS'};                         % Penalties
control.p_weights=10e-7;                          % Penalty weights
control.method='goodwin';                         % Optimisation method
control.max_iter=100;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation mode

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Gaussian pulse pair as the initial guess
pulse=[omega01_t; omega12_t];

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

end

