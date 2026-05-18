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

% Formalism and basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys);
spin_system=basis(spin_system,bas);

% Hamiltonian parameters
% Scaled Hamiltonian parameters - MHz & us
deltas=2*pi*0.5;                      % TBD - READING PAPER FOR SPECIFICS
gammas=2*pi*0;
Omega_0=2*pi*17.7;
alphas=2*pi*-295.1;


% Build the intial control pulses (FROG)
% Scaled Hamiltonian parameters - MHz & us
N=224;   % No. of Time Steps
t_g=112e-3;
t=linspace(0, t_g, N);
Fa=[-0.6137 -0.0247 0.0742 0.0507 0.0149];
Fb=[-0.0106 0.0334 0.0579 0.0140 -0.0416];

%Fa=[-0.0643    0.0718   -0.0754    0.0486    0.0234    0.0043   -0.0155    0.0105   -0.0029   -0.0683];
%Fb=[ 0.1157    0.0494   -0.0031    0.0134    0.0905    0.0027    0.0248   -0.0234   -0.0620    0.0279];

Omega_x=zeros(1,N);
Omega_y=zeros(1,N);

for n=1:length(Fa)
    Omega_x=Omega_x+Fa(n)*sin((((2*n)-1)*pi*t)/t_g);
    Omega_y=Omega_y+Fb(n)*sin((2*n*pi*t)/t_g);
end

Omega_x=Omega_0*Omega_x;
Omega_y=Omega_0*Omega_y;

% Build the Hamiltonian
H=sparse([0 0 0; 0 deltas 0; 0 0 (2*deltas+alphas)]);

% Build Control Operators
scale=(1+(gammas/Omega_0));
H_x=0.5*scale*sparse([0 1 0; 1 0 sqrt(2); 0 sqrt(2) 0]/2);
H_y=0.5*scale*sparse([0 -1i 0; 1i 0 -1i*sqrt(2); 0 1i*sqrt(2) 0]/2);

% Build source and destination states - TARGET GATE = X_pi/2
rho_init = state(spin_system,'BL1',1);   % |0>
rho_targ = (state(spin_system,'BL1',1)-1i*state(spin_system,'BL2',1))/sqrt(2);   % (|0>-i|1>)/sqrt(2)

rho_init=rho_init/norm(rho_init,'fro');
rho_targ=(rho_targ*sqrt(2))/norm(rho_targ,'fro');

% Unit fidelity is Sorensen bound
%rho_targ=rho_targ/sorensen(rho_init,rho_targ);

% Define control parameters
control.drifts={{hilb2liouv(H,'comm')}};                  % Drift
control.operators={hilb2liouv(H_x,'comm'),hilb2liouv(H_y,'comm')};     % Control
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

% Initial guess - random
pulse=[Omega_x; Omega_y];

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

% Waveform basis coefficient guess - random
%guess=randn(2,10)/20;

% Run LBFGS GRAPE pulse optimisation - Waveform basis case
%basis_coeffs=fmaxnewton(spin_system,@grape_xy,guess);

end