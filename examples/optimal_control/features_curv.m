% A transfer of coherence from longitudinal magnetization into a two-spin
% singlet state with a distribution of B1 powers. An ensemble of ten spin
% systems with different power levels is simultaneously driven to optimal
% fidelity, which in this case is 1/sqrt(2) = 0.7071
%
% Curvilinear GRAPE interface is used - the user specifies the definition
% of the curvilinear coordinates and the Jacobian. In this case, the coor-
% dinates are phase-amplitude.
%
% Calculation time: minutes.
% 
% i.kuprov@soton.ac.uk

function features_curv()

% Magnetic field
sys.magnet=14.1;

% Isotopes
sys.isotopes={'13C','13C'};

% Interactions
inter.zeeman.scalar={0.00 0.25};
inter.coupling.scalar=cell(2);
inter.coupling.scalar{1,2}=60.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,'Lz','all');
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=singlet(spin_system,1,2);
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get the control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={Lx,Ly};                        % Controls
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pulse_dt=1e-3*ones(1,50);                 % Pulse duration
control.pwr_levels=2*pi*100*linspace(0.6,1.4,11); % Power levels
control.penalties={'SNS'};                        % Penalty
control.p_weights=100;                            % Penalty weight
control.method='lbfgs';                           % Optimisation method
control.max_iter=100;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation

% Diagnostic output
control.plotting={'correlation_order','coherence_order',...
                  'xy_controls','robustness','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Curvilinear coordinates
u2x=@(u)[u(1)*cos(u(2))
         u(1)*sin(u(2))];

% Curvilinear coordinate Jacobian
dx_du=@(u)[      cos(u(2))      sin(u(2))
           -u(1)*sin(u(2)) u(1)*cos(u(2))];

% Initial guess
guess_u=[ones(1,50); (pi/2)*randn(1,50)];
       
% Run LBFGS GRAPE pulse optimisation
curv_profile=fminnewton(spin_system,@(x,y)grape_curv(x,u2x,dx_du,y),guess_u);

% Get Cartesian components of the pulse
cart_profile=zeros(size(curv_profile));
for n=1:size(curv_profile,2)
    cart_profile(:,n)=u2x(curv_profile(:,n));
end
cart_profile=mean(control.pwr_levels)*cart_profile;
CLx=cart_profile(1,:); CLy=cart_profile(2,:);

% Simulate the optimised pulse
rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{CLx,CLy},...
                    control.pulse_dt,rho_init,'expv-pwc');

% Filter out the singlet state
rho=coherence(spin_system,rho,{{'13C',0}});
rho=correlation(spin_system,rho,2,'13C');

% Run state diagnostics
stateinfo(spin_system,rho,5);

end

