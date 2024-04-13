% A transfer of coherence from a two-proton singlet state to a nearby
% carbon in a setting typically encountered in parahydrogenation expe-
% riments. LBFGS-GRAPE algorithm is used as described in
%
%            http://dx.doi.org/10.1016/j.jmr.2011.07.023
%
% Terminal fidelity is 50% in this case.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function state_transfer_s2m()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'1H','1H','13C','13C'};

% Interactions
inter.zeeman.scalar={1.5, 2.0, 30.0, 40.0};
inter.coupling.scalar=cell(4);
inter.coupling.scalar{1,2}=7.0;
inter.coupling.scalar{1,3}=150;
inter.coupling.scalar{2,4}=150;
inter.coupling.scalar{3,4}=50;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none'; 

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=singlet(spin_system,1,2);
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=state(spin_system,{'Lz'},{4});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get the control operators
LxH=operator(spin_system,'Lx','1H');
LyH=operator(spin_system,'Ly','1H');
LxC=operator(spin_system,'Lx','13C');
LyC=operator(spin_system,'Ly','13C');

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Transmitter offsets
parameters.spins={'1H','13C'};
parameters.offset=[1050 5285];
H=frqoffset(spin_system,H,parameters);

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={LxH,LyH,LxC,LyC};              % Controls
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*[460 480 500 520 540];    % Pulse power
control.pulse_dt=1.5e-4*ones(1,100);              % Slice durations
control.penalties={'NS','SNS'};                   % Penalties
control.p_weights=[0.1 10];                       % Penalty weights
control.method='lbfgs';                           % Optimisation method
control.max_iter=100;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation

% Plots during optimisation
control.plotting={'correlation_order','coherence_order',...
                  'xy_controls','local_each_spin',...
                  'spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=randn(4,100)/10;

% Run the optimisation, get normalised pulse
pulse=fminnewton(spin_system,@grape_xy,guess);

% Apply power level scaling
pulse=mean(control.pwr_levels)*pulse;
pulse=mat2cell(pulse,[1 1 1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,H,control.operators,pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');
fidelity=real(rho_targ'*rho);
report(spin_system,['Re[<target|rho(T)>] = ' num2str(fidelity)]);

end

