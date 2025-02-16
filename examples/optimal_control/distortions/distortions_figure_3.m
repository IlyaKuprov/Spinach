% Figure 3 from the paper by Rasulov and Kuprov:
%
%       https://arxiv.org/abs/2502.02198
%
% u.rasulov@soton.ac.uk
% ilya.kuprov@weizmann.ac.il

function distortions_figure_3()

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
control.pulse_dt=1e-6*ones(1,45);                  % Pulse interval grid
control.pwr_levels=2*pi*linspace(50e3,70e3,10);    % Power levels, per channel
control.method='lbfgs';                            % Optimisation method
control.penalties={'NS','SNS'};                    % Penalty types
control.p_weights=[0.01 10.0];                     % Penalty weights
control.max_iter=200;                              % Termination condition
control.parallel='ensemble';                       % Parallelisation

% Last five slices are dead time
control.freeze=zeros(2,45);
control.freeze(:,41:45)=1;

% Plotting options
control.plotting={'xy_controls','robustness','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=(1/4)*ones(2,45);
guess(:,41:45)=0;

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
figure(2); scale_figure([2.00 1.25]); subplot(2,2,3);
plot_1d(spin_system,real(spectrum),parameters);
kylabel('intensity, a.u.'); ylim([-50 210]);
klegend('GRAPE (no distortions)'); drawnow;

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
figure(2); subplot(2,2,1); 
plot_1d(spin_system,real(spectrum),parameters);
kylabel('intensity, a.u.'); ylim([-50 210]);
klegend('best square pulse'); drawnow;


% Apply RLC distortion to GRAPE pulse
omega=-sys.magnet*spin('13C'); Q=1000;
p1=exp(-abs(omega)*control.pulse_dt(1)/(2*Q));
p2=exp(-abs(omega)*control.pulse_dt(1)/(2*Q));
xy_profile=spf(xy_profile,p1);
xy_profile=spf(xy_profile,p2);
CLx=xy_profile(1,:); CLy=xy_profile(2,:);

% Simulate the distorted pulse
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
figure(2); subplot(2,2,2);
plot_1d(spin_system,real(spectrum),parameters);
kylabel('intensity, a.u.'); ylim([-50 210]);
klegend('GRAPE (RLC distorted)'); drawnow;

% Add RLC distortion to the GRAPE process
control.distortion={@(w)spf(w,p1),@(w)spf(w,p2)};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess
guess=(1/4)*ones(2,45);
guess(:,41:45)=0;

% GRAPE optimisation
xy_profile=fminnewton(spin_system,@grape_xy,guess);

% Extract the waveform
xy_profile=mean(control.pwr_levels)*xy_profile;

% Apply RLC distortion to GRAPE pulse
xy_profile=spf(xy_profile,p1);
xy_profile=spf(xy_profile,p2);
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
figure(2); subplot(2,2,4);
plot_1d(spin_system,real(spectrum),parameters);
kylabel('intensity, a.u.'); ylim([-50 210]);
klegend('RAW-GRAPE (RLC distorted)'); drawnow;

end

