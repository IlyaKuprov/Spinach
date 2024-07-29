% GRAPE Hessian internal consistency test: Newton 
% against Goodwin algorithm.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function dirdiff_5_rect()

% Set the magnetic field
sys.magnet=28.18;

% Put 100 non-interacting spins at equal intervals 
% within the [-100,+100] ppm chemical shift range 
n_spins=100; sys.isotopes=cell(n_spins,1);
for n=1:n_spins
    sys.isotopes{n}='13C';
end
inter.zeeman.scalar=num2cell(linspace(-100,100,n_spins));

% Select a basis set - IK-2 keeps complete basis on each 
% spin in this case, but ignores multi-spin orders
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.space_level=1;
bas.connectivity='scalar_couplings';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up spin states
Sp=state(spin_system,'L+','13C');
Sm=state(spin_system,'L-','13C');
Sz=state(spin_system,'Lz','13C');
Sx=(Sp+Sm)/2; Sy=(Sp-Sm)/2i;
Sx=Sx/norm(full(Sx),2);
Sy=Sy/norm(full(Sy),2);
Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lp=operator(spin_system,'L+','13C');
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={Lx,Ly};                      % Controls
control.rho_init={ Sx Sy Sz};                   % Starting states
control.rho_targ={-Sz Sy Sx};                   % Target states
control.pwr_levels=2*pi*linspace(50e3,70e3,10); % Power levels
control.max_iter=1000;                          % Termination condition
control.plotting={};                            % Plotting options
control.integrator='rectangle';                 % Integrator
control.pulse_dt=12.8e-6*ones(1,5);             % Time interval grid

% Pick initial guess, phase-modulated GRAPE
control.amplitudes=ones(1,5); guess=randn(1,5)/3;

% Get Newton Hessian
control.method='newton';
spin_system=optimcon(spin_system,control);
[~,~,~,newton_hess_ph]=grape_phase(guess,spin_system);

% Get Goodwin Hessian
control.method='goodwin';
spin_system=optimcon(spin_system,control);
[~,~,~,goodwin_hess_ph]=grape_phase(guess,spin_system);

% Pick initial guess, XY-modulated GRAPE
guess=randn(2,5)/3;

% Get Newton Hessian
control.method='newton';
spin_system=optimcon(spin_system,control);
[~,~,~,newton_hess_xy]=grape_xy(guess,spin_system);

% Get Goodwin Hessian
control.method='goodwin';
spin_system=optimcon(spin_system,control);
[~,~,~,goodwin_hess_xy]=grape_xy(guess,spin_system);

% Run the comparisons
if norm(newton_hess_ph(:)-goodwin_hess_ph(:),1)>1e-6*norm(newton_hess_ph(:),1)
    error('Phase Hessian internal consistency test failed.');
else
    disp('Phase Hessian internal consistency test passed.');
end
if norm(newton_hess_xy(:)-goodwin_hess_xy(:),1)>1e-6*norm(newton_hess_xy(:),1)
    error('Cartesian Hessian internal consistency test failed.');
else
    disp('Cartesian internal consistency test passed.');
end

end

