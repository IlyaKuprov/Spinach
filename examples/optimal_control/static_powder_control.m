% Optimal control optimisation for a pulse that is designed 
% to set deuterium magnetisation in a -CD3 group of alanine
% up for perfect rephasing 100 microseconds after the pulse
% is finished.
%
% The system is a powder (100 orientations) with a B1 dist-
% ribution (from 46 to 54 kHz per channel) and transmitter
% offset error within 1 kHz of the chemical shift.
%
% Goodwin's very efficient version of the GRAPE Hessian al-
% gorithm is used because propagator dimensions are small;
% it yields a sophisticated kind of spin echo.
%
% Calculation time: minutes
%
% ilya.kuprov@weizmann.ac.il
% u.rasulov@soton.ac.uk
% marina.carravetta@soton.ac.uk

function static_powder_control()

% 600 MHz magnet
sys.magnet=14.1; 

% Isotopes
sys.isotopes={'2H'};

% Alanine CD3 NQI parameters
inter.coupling.matrix{1,1}=anas2mat(0,40e3,0,0,0,0);

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none'; 

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Powder context parameters
parameters.spins={'2H'};
parameters.decouple={};
parameters.offset=0;
parameters.grid='rep_2ang_100pts_sph';
parameters.rframes={{'2H',2}};
parameters.verbose=0;

% Drift Liouvillians for the ensemble
control.drifts=drifts(spin_system,@powder,parameters,'labframe');

% Set up and normalise the initial state
rho_init=state(spin_system,'Lz','2H');
rho_init=rho_init/norm(full(rho_init),2);

% Set up and normalise the target state
rho_targ=state(spin_system,'Lx','2H');
rho_targ=rho_targ/norm(full(rho_targ),2);

% Get deuterium control operators
Dx=operator(spin_system,'Lx','2H');
Dy=operator(spin_system,'Ly','2H');

% Get deuterium offset operator
Dz=operator(spin_system,'Ly','2H');

% Define control parameters
control.operators={Dx,Dy};                     % Controls
control.rho_init={rho_init};                   % Starting state
control.rho_targ={rho_targ};                   % Destination state
control.pwr_levels=2*pi*[46 48 50 52 54]*1e3;  % Power distribution
control.off_ops={Dz};                          % Offset operator
control.offsets={linspace(-1e3,1e3,5)};        % Offset distribution
control.pulse_dt=2e-6*ones(1,100);             % Slice duration grid
control.method='goodwin';                      % Optimisation method
control.max_iter=100;                          % Maximum iterations
control.dead_time=100e-6;                      % Dead time
control.penalties={'NS','SNS'};                % Penalties
control.p_weights=[0.1 10];                    % Penalty weights
control.parallel='ensemble';                   % Parallelisation

% Plotting options
control.plotting={'xy_controls','robustness','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Waveform guess
guess=randn(2,100)/10;

% Run the optimisation
pulse_profile=fminnewton(spin_system,@grape_xy,guess);

% Get Cartesian components of the pulse
CLx=mean(control.pwr_levels)*pulse_profile(1,:);
CLy=mean(control.pwr_levels)*pulse_profile(2,:);

% Run a test calculation
fid_optim=zeros([500 1],'like',1i); 
fid_ideal=zeros([500 1],'like',1i);
parfor n=1:numel(control.drifts)
    
    % Run the pulse
    rho=shaped_pulse_xy(spin_system,control.drifts{n}{1},{Dx,Dy},{CLx,CLy},...
                        control.pulse_dt,rho_init,'expv-pwc'); %#ok<PFBNS>
    
    % Run the evolution
    fid_optim=fid_optim+evolution(spin_system,control.drifts{n}{1},...
                                  rho_targ,rho,0.5e-6,499,'observable');
    fid_ideal=fid_ideal+evolution(spin_system,control.drifts{n}{1},...
                                  rho_targ,rho_targ,0.5e-6,499,'observable');
    
end

% Plot the echo
figure(); subplot(1,2,1);
time_axis=linspace(0,250,500)';
plot(time_axis,real(fid_optim)); 
kxlabel('time after pulse, $\mu$s');
grid on; xlim('tight');

% The OC echo and the ideal FID
fid_optim=fid_optim(201:end);
fid_ideal=fid_ideal(1:300);
fid_optim=apodisation(spin_system,fid_optim,{{'exp',5}});
fid_ideal=apodisation(spin_system,fid_ideal,{{'exp',5}});
spec_optim=real(fftshift(fft(fid_optim)));
spec_ideal=real(fftshift(fft(fid_ideal)));

% Plot Fourier transform comparison
freq_axis=linspace(-1/0.5e-6,1/0.5e-6,300)'/2000;
subplot(1,2,2); grid on; axis tight; hold on;
plot(freq_axis,[spec_optim/max(spec_optim) ...
                spec_ideal/max(spec_ideal)]);
klegend({'half-echo fft','ideal fid fft'});
xlim([-200 200]); kxlabel('frequency, kHz');

end

