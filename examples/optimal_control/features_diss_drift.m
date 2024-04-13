% Optimal control optimisation of a pulse performing magnetisa-
% tion transfer from H(N) to C(O) in a typical protein backbone
% spin system (literature data for shifts and couplings) with a 
% range of pulse powers emulating B1 inhomogeneity.
%
% The dynamics includes dissipative terms in the drift generator:
% C(O) and N(H) are set to have rapid transverse relaxation. Four-
% spin correlation approximation is used, wherein five-spin and 
% higher correlations are dropped from the basis set.
%
% Calculation time: hours.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function features_diss_drift()

% Magnetic field
sys.magnet=9.4;

% Spin system
sys.isotopes={'15N','1H','13C','13C','13C','15N'};
sys.labels={'N_(n)','H','CA','CB','C','N_(n+1)'};

% Textbook chemical shifts, ppm
inter.zeeman.scalar={119.79, 8.03, 57.32, 27.71, 177.25, 115.55};

% Scalar couplings, Hz (literature values)
inter.coupling.scalar=cell(6);
inter.coupling.scalar{1,3}=-11; 
inter.coupling.scalar{2,3}=140;
inter.coupling.scalar{3,4}=35;
inter.coupling.scalar{3,5}=55;
inter.coupling.scalar{3,6}=7;
inter.coupling.scalar{5,6}=-15;

% Relaxation theory
inter.relaxation={'t1_t2'};
inter.r1_rates={ 1.0 1.0 1.0 1.0   1.0  1.0};
inter.r2_rates={50.0 1.0 1.0 1.0 100.0 50.0};
inter.equilibrium='zero';
inter.rlx_keep='diagonal';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=4;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,{'Lz'},{2});
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=state(spin_system,{'Lz'},{5});
rho_targ=rho_targ/norm(full(rho_targ),2);

% Control operators
LxH=operator(spin_system,'Lx','1H');
LxC=operator(spin_system,'Lx','13C');
LxN=operator(spin_system,'Lx','15N');
LyH=operator(spin_system,'Ly','1H');
LyC=operator(spin_system,'Ly','13C');
LyN=operator(spin_system,'Ly','15N');

% Dissipative drift Liouvillian
spin_system=assume(spin_system,'nmr');
L=hamiltonian(spin_system)+1i*relaxation(spin_system);

% Put transmitters in the right place
parameters.spins={'1H','13C','15N'};
parameters.offset=[3214 10000 -4800];
L=frqoffset(spin_system,L,parameters);

% Define control parameters
control.drifts={{L}};                             % Drift
control.operators={LxH,LyH,LxC,LyC,LxN,LyN};      % Controls
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*linspace(0.9e3,1.1e3,5);  % Pulse power, rad/s
control.pulse_dt=5e-5*ones(1,400);                % Slice durations
control.penalties={'NS','SNS'};                   % Penalties
control.p_weights=[0.1 10];                       % Penalty weights
control.method='lbfgs';                           % Optimisation method
control.max_iter=500;                             % Termination tolerance
control.parallel='ensemble';                      % Parallelisation

% Control trajectory analysis plots
control.plotting={'correlation_order','local_each_spin',...
                  'amp_controls','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Start with pi/2 on 1H;
% end with pi/2 on 13C
guess=randn(6,400)/10; 
guess(2,1:20)=0.25;
guess(4,(end-20):end)=0.25;

% Run the optimisation
pulse=fminnewton(spin_system,@grape_xy,guess);

% Apply power level scaling
pulse=mean(control.pwr_levels)*pulse;
pulse=mat2cell(pulse,[1 1 1 1 1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,L,control.operators,pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');
fidelity=real(rho_targ'*rho);
report(spin_system,['Re[<target|rho(T)>] = ' num2str(fidelity)]);

end

