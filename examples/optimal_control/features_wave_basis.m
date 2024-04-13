% Optimal control optimisation of a pulse performing magnetisation
% transfer from proton 2 to carbon 5 in a typical protein backbone
% spin system (literature data for chemical shifts, with suitable 
% transmitter offsets, and couplings).
%
% The waveform is optimised with the standard LBFGS-GRAPE algorithm
% in a sine wave basis set with a penalty on the waveform exceeding
% a defined power threshold.
%
% Spin correlation order populations and the waveform are plotted
% at every iteration, the fidelity is printed to the console.
%
% Note: optimal control problems using waveform basis sets (as op-
%       posed to varying the waveform directly point-by-point) ha-
%       ve slow convergence rates -- about 1000 iterations are re-
%       quired to achieve a good fidelity in this case.
%
% http://dx.doi.org/10.1016/j.jmr.2011.07.023
%
% Calculation time: hours.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function features_wave_basis()

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

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Put transmitters in the right place
parameters.spins={'1H','13C','15N'};
parameters.offset=[3214 10000 -4800];
H=frqoffset(spin_system,H,parameters);

% Define control parameters
control.drifts={{H}};	                        % Drift
control.operators={LxH,LyH,LxC,LyC,LxN,LyN};    % Controls
control.rho_init={rho_init};                    % Starting state
control.rho_targ={rho_targ};                    % Destination state
control.pwr_levels=2*pi*linspace(10e3,12e3,5);  % Pulse power, rad/s
control.pulse_dt=1e-4*ones(1,200);              % Pulse duration
control.basis=wave_basis('sine_waves',30,200)'; % Basis set
control.penalties={'SNS'};                      % Penalty
control.p_weights=100;                          % Penalty weight
control.method='lbfgs';                         % Optimisation method
control.max_iter=1000;                          % Termination tolerance
control.parallel='ensemble';                    % Parallelisation
 
% Control trajectory analysis plots
control.plotting={'correlation_order','local_each_spin',...
                  'xy_controls','robustness','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=randn(6,30);    

% Run the optimisation
basis_coeffs=fminnewton(spin_system,@grape_xy,guess);

% Reassemble time-domain control sequence
pulse=mean(control.pwr_levels)*basis_coeffs*control.basis;
pulse=mat2cell(pulse,[1 1 1 1 1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,H,control.operators,pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');
fidelity=real(rho_targ'*rho);
report(spin_system,['Re[<target|rho(T)>] = ' num2str(fidelity)]);

end

