% Directional derivative test for the GRAPE with Distortions
%
% i.kuprov@soton.ac.uk
% u.rasulov@soton.ac.uk

function dirdiff_7_distortions()

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

% Get the control operator
Lp=operator(spin_system,'L+','13C');
Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i;

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                           % Drift
control.operators={Lx Ly};                      % Controls
control.rho_init={ Sx Sy Sz};                   % Starting states
control.rho_targ={-Sz Sy Sx};                   % Target states
control.pwr_levels=2*pi*linspace(50e3,70e3,10); % Power levels
control.method='lbfgs';                         % Optimisation method
control.max_iter=1000;                          % Termination condition
control.plotting={};                            % Plotting options
control.integrator='rectangle';                 % Integrator

n_steps=10;
h=sqrt(eps('double'));

% Set the interval grid
control.pulse_dt=12.8e-6*ones(1,n_steps);

% ------------ TEST NO DISTORTIONS RESPONSE FUNCTION -------------------
responses_applied={@(distorted_waveform)no_dist(distorted_waveform)};
control.distortion = responses_applied;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Random guess and finite diff increment
guess=ones(2,n_steps)/3; 

% Call GRAPE and request analytical gradient
[~,~,grad_anl]=grape_xy(guess,spin_system);
grad_anl=squeeze(grad_anl(:,:,1));

% Left waveform edge
wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
wave_back=guess; wave_back(1)=wave_back(1)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6
    disp('left edge test passed');
else
    error('left edge test failed');
end

% Right waveform edge
wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
wave_back=guess; wave_back(end)=wave_back(end)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6
    disp('right edge test passed');
else
    error('right edge test failed');
end

% Waveform midpoint
wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
wave_back=guess; wave_back(3)=wave_back(3)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6
    disp('midpoint test passed');
else
    error('midpoint test failed');
end

% ------------ TEST RECURSIVE FILTER RESPONSE FUNCTION -------------------
% Set Recursive Filter Parameters
 a = zeros(3,1); b = zeros(2,1);
 a(1) = 0.05; a(2)=0.2; a(3)=0; b(1)=0.95; b(2)=0.6;

responses_applied={@(distorted_waveform)recursive_filter(distorted_waveform,a,b)};
control.distortion = responses_applied;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Random guess and finite diff increment
guess=ones(2,n_steps)/3; 

% Call GRAPE and request analytical gradient
[~,~,grad_anl]=grape_xy(guess,spin_system);
grad_anl=squeeze(grad_anl(:,:,1));

% Left waveform edge
wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
wave_back=guess; wave_back(1)=wave_back(1)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6
    disp('left edge test passed');
else
    error('left edge test failed');
end

% Right waveform edge
wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
wave_back=guess; wave_back(end)=wave_back(end)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6
    disp('right edge test passed');
else
    error('right edge test failed');
end

% Waveform midpoint
wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
wave_back=guess; wave_back(3)=wave_back(3)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6
    disp('midpoint test passed');
else
    error('midpoint test failed');
end


% ------------ TEST 3 SEQUENTIAL RESPONSE FUNCTIONS -------------------
% Set Recursive Filter Parameters
 a = zeros(3,1); b = zeros(2,1);
 a(1) = 0.05; a(2)=0.2; a(3)=0; b(1)=0.95; b(2)=0.6;

% Set Amplifier Saturation Parameters
saturation_level=1; multiplier=mean(2 * pi * linspace(50e3, 70e3, 10));
saturation=saturation_level*multiplier;

responses_applied={@(distorted_waveform)amp_comp(distorted_waveform,saturation) ...
,@(distorted_waveform)no_dist(distorted_waveform),@(distorted_waveform)recursive_filter(distorted_waveform,a,b)};

control.distortion = responses_applied;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Random guess and finite diff increment
guess=ones(2,n_steps)/3; 

% Call GRAPE and request analytical gradient
[~,~,grad_anl]=grape_xy(guess,spin_system);
grad_anl=squeeze(grad_anl(:,:,1));

% Left waveform edge
wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
wave_back=guess; wave_back(1)=wave_back(1)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6
    disp('left edge test passed');
else
    error('left edge test failed');
end

% Right waveform edge
wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
wave_back=guess; wave_back(end)=wave_back(end)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6
    disp('right edge test passed');
else
    error('right edge test failed');
end

% Waveform midpoint
wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
wave_back=guess; wave_back(3)=wave_back(3)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6
    disp('midpoint test passed');
else
    error('midpoint test failed');
end

% ------------ TEST ENSEMBLE OF RESPONSE FUNCTIONS -------------------

% Set Amplifier Saturation Parameters
saturation_level=1; multiplier=mean(2 * pi * linspace(50e3, 70e3, 10));

responses_applied={@(distorted_waveform)amp_comp(distorted_waveform,saturation),@(distorted_waveform)no_dist(distorted_waveform); ...
@(distorted_waveform)amp_comp(distorted_waveform,(saturation_level+1)*multiplier),@(distorted_waveform)no_dist(distorted_waveform);...
@(distorted_waveform)amp_comp(distorted_waveform,(saturation_level+2)*multiplier),@(distorted_waveform)no_dist(distorted_waveform);...
@(distorted_waveform)amp_comp(distorted_waveform,(saturation_level+3)*multiplier),@(distorted_waveform)no_dist(distorted_waveform)};

control.distortion = responses_applied;

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Random guess and finite diff increment
guess=ones(2,n_steps)/3; 

% Call GRAPE and request analytical gradient
[~,~,grad_anl]=grape_xy(guess,spin_system);
grad_anl=squeeze(grad_anl(:,:,1));

% Left waveform edge
wave_forw=guess; wave_forw(1)=wave_forw(1)+h;
wave_back=guess; wave_back(1)=wave_back(1)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6
    disp('left edge test passed');
else
    error('left edge test failed');
end

% Right waveform edge
wave_forw=guess; wave_forw(end)=wave_forw(end)+h;
wave_back=guess; wave_back(end)=wave_back(end)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6
    disp('right edge test passed');
else
    error('right edge test failed');
end

% Waveform midpoint
wave_forw=guess; wave_forw(3)=wave_forw(3)+h;
wave_back=guess; wave_back(3)=wave_back(3)-h;
[~,fid_forw]=grape_xy(wave_forw,spin_system);
[~,fid_back]=grape_xy(wave_back,spin_system);
grad_num=(fid_forw(1)-fid_back(1))/(2*h);
if abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6
    disp('midpoint test passed');
else
    error('midpoint test failed');
end

end




