% Basic implementation of a Frequency Robust Gate 
% (FROG) for a single transmon, based on: 
%
%   https://doi.org/10.48550/arXiv.2511.22580
%
% Ensemble GRAPE optimisation with a distribution
% of control powers and excess amplitude penalty.
%
% Calculation time: minutes
%
% c.musselwhite@soton.ac.uk

function transmon_frog()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3'};

% Rotating frame transmon parameters
inter.modes.frqs={0.5e6};
inter.modes.anharms={-295.1e6};

% Formalism and basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Drift Hamiltonian from the declared interactions
H=hamiltonian(assume(spin_system,'cavity'));

% Build the intial control pulses (FROG)
nsteps=224; t_g=112e-9; t=linspace(0,t_g,nsteps);
Fa=[-0.6137 -0.0247 0.0742 0.0507 0.0149];
Fb=[-0.0106 0.0334 0.0579 0.0140 -0.0416];
omega_x=zeros(1,nsteps); omega_y=zeros(1,nsteps);
for n=1:length(Fa)
    omega_x=omega_x+Fa(n)*sin((((2*n)-1)*pi*t)/t_g);
    omega_y=omega_y+Fb(n)*sin((2*n*pi*t)/t_g);
end

% Quadrature control operators
Cr=operator(spin_system,'C',1); 
An=operator(spin_system,'A',1);
Cx=(Cr+An)/2; Cy=1i*(Cr-An)/2;

% Source and destination states
psi_targ=[1; -1i; 0]/sqrt(2);
rho_init=state(spin_system,'BL1',1);
rho_targ=hilb2liouv(psi_targ*psi_targ','statevec');
rho_init=rho_init/norm(rho_init,2);
rho_targ=rho_targ/norm(rho_targ,2);

% Define control parameters
control.isotopes={'T3'};                             % Isotopes
control.channels=[1;1];                              % Channel map
control.drifts={{H}};                                % Drift
control.operators={Cx,Cy};                           % Control
control.rho_init={rho_init};                         % Starting state
control.rho_targ={rho_targ};                         % Destination state
control.pwr_levels=2*pi*[15e6 16e6 17e6 18e6 19e6];  % Pulse power ensemble 
control.pulse_dt=(t_g/nsteps)*ones(1,nsteps);        % Slice durations
control.penalties={'SNSA'};                          % Penalties
control.p_weights=0.01;                              % Penalty weights
control.method='goodwin';                            % Optimisation method
control.max_iter=50;                                 % Termination condition

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
pulse=[omega_x; omega_y];

% Run the optimisation
fmaxnewton(spin_system,@grape_xy,pulse);

end

