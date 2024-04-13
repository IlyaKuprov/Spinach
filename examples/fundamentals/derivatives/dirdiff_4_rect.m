% Directional derivative test for the phase-modulated GRAPE
% module, rectangles integrator.
%
% i.kuprov@soton.ac.uk
% uluk.rasulov@soton.ac.uk

function dirdiff_4_rect()

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
control.method='lbfgs';                         % Optimisation method
control.max_iter=1000;                          % Termination condition
control.plotting={};                            % Plotting options
control.integrator='rectangle';                 % Integrator

% Set the interval grid
control.pulse_dt=12.8e-6*ones(1,5);
control.amplitudes=ones(1,5); 

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Random phases and finite diff increment
guess=randn(1,5)/3; h=sqrt(eps('double'));

% Call GRAPE and request analytical gradient
[~,~,grad_anl]=grape_phase(guess,spin_system);
grad_anl=squeeze(grad_anl(:,:,1));

% Left waveform edge
wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
wave_back=guess; wave_back(1)=wave_back(1)-h;
[~,fid_forw]=grape_phase(wave_forw,spin_system);
[~,fid_back]=grape_phase(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6
    disp('left edge test passed');
else
    error('left edge test failed');
end

% Right waveform edge
wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
wave_back=guess; wave_back(end)=wave_back(end)-h;
[~,fid_forw]=grape_phase(wave_forw,spin_system);
[~,fid_back]=grape_phase(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6
    disp('right edge test passed');
else
    error('right edge test failed');
end

% Waveform midpoint
wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
wave_back=guess; wave_back(3)=wave_back(3)-h;
[~,fid_forw]=grape_phase(wave_forw,spin_system);
[~,fid_back]=grape_phase(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6
    disp('midpoint test passed');
else
    error('midpoint test failed');
end

end

