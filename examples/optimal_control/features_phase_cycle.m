% Optimal control pulse optimisation for state-to-state transfer across a
% scalar coupling in a hydrofluorocarbon fragment spin system. The start-
% ing state is Z-magnetisation on 1H, the destination state is quadrature
% transverse magnetisation on 19F. 
%
% A phase cycle is specified: a flip in the phase of the fluorine channel
% must produce the corresponding flip in the phase of the resulting magne-
% tisation on 19F.
%
% Calculation time: minutes.
%
% david.goodwin@inano.au.dk
% ilya.kuprov@weizmann.ac.il

function features_phase_cycle()

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
rho_targ=(state(spin_system,{'L+'},{3})+...
          state(spin_system,{'L-'},{3}));
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
control.pulse_dt=2e-4*ones(1,50);                 % Slice durations
control.penalties={'SNS'};                        % Penalty
control.p_weights=100;                            % Penalty weight
control.method='lbfgs';                           % Optimisation method
control.max_iter=100;                             % Termination condition
control.phase_cycle=[0  0  0  0  0;
                     0  0  0  pi pi];             % Phase cycle (only 0, pi supported so far)
control.parallel='ensemble';                      % Parallelisation

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

% Run test simulations
for n=1:size(control.phase_cycle,1)

    % Apply the phase to the pulse
    f_channel=pulse(5:6,:); phi=control.phase_cycle(n,4);
    f_channel=[cos(phi), -sin(phi); sin(phi), cos(phi)]*f_channel;
    phased_pulse=pulse; phased_pulse(5:6,:)=f_channel;

    % Run a test simulation using the optimal pulse
    report(spin_system,['Test simulation, phase cycle step ' int2str(n) ':']);
    phased_pulse=mat2cell(phased_pulse,[1 1 1 1 1 1]);
    rho=shaped_pulse_xy(spin_system,H,control.operators,phased_pulse,...
                                      control.pulse_dt,rho_init,'expv-pwc');

    % Run state diagnoistics
    stateinfo(spin_system,rho,10);

end

end

