% Basic implementation of a Frequency Robust Gate (FROG)
% for a single Transmon qubit, based on: (https://doi.org/10.48550/arXiv.2511.22580)
% by S.Glaser et al.
%
% c.musselwhite@soton.ac.uk
%

function two_transmons_robust()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3'};

% Rotating frame transmon parameters
% Scaled Hamiltonian parameters - MHz & us
inter.modes.frqs={0.5};
inter.modes.anharms={-295.1};

% Formalism and basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Build the intial control pulses (FROG)
% Scaled Hamiltonian parameters - MHz & us
N=224;   % No. of Time Steps
t_g=112e-3;
t=linspace(0, t_g, N);
Fa=[-0.6137 -0.0247 0.0742 0.0507 0.0149];
Fb=[-0.0106 0.0334 0.0579 0.0140 -0.0416];

Omega_x=zeros(1,N);
Omega_y=zeros(1,N);

for n=1:length(Fa)
    Omega_x=Omega_x+Fa(n)*sin((((2*n)-1)*pi*t)/t_g);
    Omega_y=Omega_y+Fb(n)*sin((2*n*pi*t)/t_g);
end

% Quadrature control operators of the transmon mode
Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1);
Cx=(Cr+An)/2; Cy=1i*(Cr-An)/2;

% Build source and destination states - TARGET GATE = X_pi/2
psi_targ=[1; -1i; 0]/sqrt(2);
rho_init=state(spin_system,'BL1',1);
rho_targ=hilb2liouv(psi_targ*psi_targ','statevec');
rho_init=rho_init/norm(rho_init,'fro');
rho_targ=rho_targ/norm(rho_targ,'fro');

% Unit fidelity is Sorensen bound
%rho_targ=rho_targ/sorensen(rho_init,rho_targ);

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={Cx,Cy};                        % Control
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
%control.basis=[wave_basis('sine_waves',5,N) wave_basis('cosine_waves',5,N)]';    % Basis set
control.pwr_levels=2*pi*[15 16 17 18 19];       % Pulse power ensemble (ADD *1e6 WHEN NOT USING SCALED PARAMETERS)
control.pulse_dt=(t_g/N)*ones(1,N);                % Slice durations
control.penalties={'NS'};                         % Penalties
control.p_weights=1e-5;                            % Penalty weights
control.method='goodwin';                           % Optimisation method
control.max_iter=100;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation mode

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess - normalised FROG waveform, scaled by the power levels
pulse=[Omega_x; Omega_y];

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

% Waveform basis coefficient guess - random
%guess=randn(2,10)/20;

% Run LBFGS GRAPE pulse optimisation - Waveform basis case
%basis_coeffs=fmaxnewton(spin_system,@grape_xy,guess);

end

