% A template file for the "magic pulse" optimisations. The term refers to
% a family of broadband NMR pulses that are tolerant to resonance offsets
% and power calibration errors:
%
%               http://dx.doi.org/10.1016/j.jmr.2005.12.010
%
% Consider a 13C 90-degree excitation pulse in a 28.18 Tesla magnet. The
% pulse must uniformly excite a bandwidth of around 200 ppm (60 kHz) and
% must be short enough for the worst-case 13C-1H J-coupling (ca. 200 Hz)
% to be negligible. The latter requirement caps the duration at 1/100*J
% = 50 us. The pulse must accomplish the following transfers: {Lz -> Lx,
% Ly -> Ly, Lx -> -Lz}. A realistically achievable nutation frequency is
% between 50 kHz and 70 kHz across the RF coil.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk
% david.goodwin@inano.au.dk

function magic_pulse_cart()

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
Sx=state(spin_system,'Lx','13C');
Sy=state(spin_system,'Ly','13C');
Sz=state(spin_system,'Lz','13C');
Sx=Sx/norm(full(Sx),2);
Sy=Sy/norm(full(Sy),2);
Sz=Sz/norm(full(Sz),2);

% Get the control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');

% Get the drift Hamiltonian
H=hamiltonian(assume(spin_system,'nmr'));

% Define control parameters
control.drifts={{H}};                              % Drift
control.operators={Lx,Ly};                         % Controls
control.rho_init={ Sx Sy Sz};                      % Starting states
control.rho_targ={-Sz Sy Sx};                      % Target states
control.pulse_dt=1e-6*ones(1,40);                  % Pulse interval grid
control.pwr_levels=2*pi*linspace(50e3,70e3,10);    % Power levels, per channel
control.method='lbfgs';                            % Optimisation method
control.penalties={'NS','SNS'};                    % Penalty types
control.p_weights=[0.01 10.0];                     % Penalty weights
control.max_iter=200;                              % Termination condition
control.parallel='ensemble';                       % Parallelisation

% Plotting options
control.plotting={'xy_controls','robustness','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=(1/4)*ones(2,40);

% GRAPE optimisation
xy_profile=fminnewton(spin_system,@grape_xy,guess);

% Extract the waveform
xy_profile=mean(control.pwr_levels)*xy_profile;
CLx=xy_profile(1,:); CLy=xy_profile(2,:);

% Simulate the optimised pulse
rho_init=state(spin_system,'Lz','13C');
rho=shaped_pulse_xy(spin_system,H,{Lx,Ly},{CLx,CLy},...
                    control.pulse_dt,rho_init,'expv-pwc');

% Set acquisition parameters
parameters.spins={'13C'};
parameters.rho0=rho;
parameters.coil=state(spin_system,'L+','13C');
parameters.decouple={};
parameters.offset=0;
parameters.sweep=70000;
parameters.npoints=2048;
parameters.zerofill=16384;
parameters.axis_units='ppm';
parameters.invert_axis=1;

% Simulate the free induction decay
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'gauss',10}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(2); subplot(2,1,2);
plot_1d(spin_system,real(spectrum),parameters);
kylabel('intensity, a.u.');

% Simulate the conventional hard pulse
parameters.rho0=state(spin_system,'Lz','13C');
parameters.pulse_frq=0;
parameters.pulse_phi=pi/2;
parameters.pulse_pwr=2*pi*60e3;
parameters.pulse_dur=4.2e-6;
parameters.pulse_rnk=3;
parameters.method='expv';
fid=liquid(spin_system,@sp_acquire,parameters,'nmr');
fid=apodisation(spin_system,fid,{{'gauss',10}});
spectrum=fftshift(fft(fid,parameters.zerofill));
subplot(2,1,1); plot_1d(spin_system,real(spectrum),parameters);
kylabel('intensity, a.u.');

end

