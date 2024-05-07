% Optimal control optimisation of a pulse performing magnetisa-
% tion transfer from H(N) to C(O) in a typical protein backbone
% spin system (literature data for shifts and couplings) with a 
% range of pulse powers emulating B1 inhomogeneity and a range
% of offsets to account for imperfect transmitter placement.
%
% The waveform is optimized with LBFGS-GRAPE algorithm with po-
% int-by-point variation and a penalty on the pulse amplitude.
%
%         http://dx.doi.org/10.1016/j.jmr.2011.07.023
%
% Calculation time: hours.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function state_transfer_pro()

% Magnetic field
sys.magnet=9.4;

% Spin system
sys.isotopes={'15N','1H','13C','13C','13C','15N'};
sys.labels={'N_{(n)}','H','CA','CB','C','N_{(n+1)}'};

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
LyH=operator(spin_system,'Ly','1H');
LxC=operator(spin_system,'Lx','13C');
LyC=operator(spin_system,'Ly','13C');
LxN=operator(spin_system,'Lx','15N');
LyN=operator(spin_system,'Ly','15N');

% Offset operators
LzH=operator(spin_system,'Lz','1H');
LzC=operator(spin_system,'Lz','13C');
LzN=operator(spin_system,'Lz','15N');

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Put transmitters in the right place
parameters.spins={'1H','13C','15N'};
parameters.offset=[3214 10000 -4800];
H=frqoffset(spin_system,H,parameters);

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={LxH,LyH,LxC,LyC,LxN,LyN};      % Controls
control.off_ops={LzH,LzC,LzN};                    % Offset operators 
control.offsets={linspace(-100,100,3),...
                 linspace(-100,100,3),...
                 linspace(-100,100,3)};           % Offset ranges, Hz
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*linspace(0.9e3,1.1e3,5);  % Pulse power, rad/s
control.pulse_dt=4e-5*ones(1,500);                % Slice durations
control.penalties={'NS'};                         % Penalty function
control.p_weights=0.01;                           % Penalty weight
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
guess=randn(6,500)/10; 
guess(2,1:20)=0.25;
guess(4,(end-20):end)=0.25;

% Run the optimisation
pulse=fminnewton(spin_system,@grape_xy,guess);

% Apply power level scaling
pulse=mean(control.pwr_levels)*pulse;
pulse=mat2cell(pulse,[1 1 1 1 1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,H,control.operators,pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');
fidelity=real(rho_targ'*rho);
report(spin_system,['Re[<target|rho(T)>] = ' num2str(fidelity)]);

end

