% GRAPE Hessian test against finite-differenced gradients.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function dirdiff_6_rect()

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
control.method='newton';                        % Optimisation method

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Random guess and finite diff increment
guess=randn(2,5)/3; h=sqrt(eps('double'));

% Call GRAPE and request analytical Hessian
[~,~,~,hess_anl]=grape_xy(guess,spin_system);
hess_anl=squeeze(hess_anl(:,:,1));

% Leftmost Hessian column
wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
wave_back=guess; wave_back(1)=wave_back(1)-h;
[~,~,grad_forw]=grape_xy(wave_forw,spin_system);
[~,~,grad_back]=grape_xy(wave_back,spin_system);
grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
hess_num=(grad_forw-grad_back)/(2*h);
if norm(hess_anl(:,1)-hess_num,1)<1e-6*norm(hess_num,1)
    disp('leftmost column test passed');
else
    error('leftmost column test failed');
end

% Rightmost Hessian column
wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
wave_back=guess; wave_back(end)=wave_back(end)-h;
[~,~,grad_forw]=grape_xy(wave_forw,spin_system);
[~,~,grad_back]=grape_xy(wave_back,spin_system);
grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
hess_num=(grad_forw-grad_back)/(2*h);
if norm(hess_anl(:,end)-hess_num,1)<1e-6*norm(hess_num,1)
    disp('rightmost column test passed');
else
    error('rightmost column test failed');
end

% Middle Hessian column
wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
wave_back=guess; wave_back(3)=wave_back(3)-h;
[~,~,grad_forw]=grape_xy(wave_forw,spin_system);
[~,~,grad_back]=grape_xy(wave_back,spin_system);
grad_forw=squeeze(grad_forw(:,:,1)); grad_forw=grad_forw(:);
grad_back=squeeze(grad_back(:,:,1)); grad_back=grad_back(:);
hess_num=(grad_forw-grad_back)/(2*h);
if norm(hess_anl(:,3)-hess_num,1)<1e-6*norm(hess_num,1)
    disp('middle column test passed');
else
    error('middle column test failed');
end

end

