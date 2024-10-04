% A transfer of coherence from longitudinal magnetization into a
% two-spin singlet state in allyl pyruvate with a distribution of
% B1 powers and transmitter offsets. XX and YY components of the 
% singlet dephase rapidly in this system, and are therefore drop-
% ped from the target state specification.
%
% Calculation time: many hours.
%
% i.kuprov@soton.ac.uk
% a.acharya@soton.ac.uk

function state_transfer_m2s()

% Get the spin system
[sys,inter]=allyl_pyruvate({'1H'});

% Kill the methyl group
sys.isotopes=sys.isotopes(1:5);
sys.labels=sys.labels(1:5);
inter.zeeman.scalar=inter.zeeman.scalar(1:5);
inter.zeeman.matrix=inter.zeeman.matrix(1:5);
inter.coordinates=inter.coordinates(1:5);
inter.coupling.scalar=inter.coupling.scalar(1:5,1:5);

% Magnetic field (500.13 MHz) 
sys.magnet=11.7464;

% Formalism and basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Set up and normalise the initial state
rho_init=state(spin_system,'Lz','all');

% Set up and normalise target state 
rho_targ=state(spin_system,{'Lz','Lz'},...
                           {idxof(sys,'Hb'),...
                            idxof(sys,'Hc')});

% Get the control operators
LxH=operator(spin_system,'Lx','1H');
LyH=operator(spin_system,'Ly','1H');

% Drift Hamiltonian 
H=hamiltonian(assume(spin_system,'nmr'));

% Transmitter offset
parameters.spins={'1H'}; parameters.offset=2670;
H=frqoffset(spin_system,H,parameters);

% Offset distribution generator
Hz=operator(spin_system,'Lz','1H');

% Define control parameters
control.drifts={{H}};                             % Drift
control.operators={LxH,LyH};                      % Controls
control.off_ops={Hz};                             % Offset operator
control.offsets={linspace(-10,10,5)};             % Offset distribution
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*[480 490 500 510 520];    % Pulse power ensemble
control.pulse_dt=1e-3*ones(1,300);                % Slice durations
control.penalties={'NS','SNS'};                   % Penalties
control.p_weights=[1 10];                         % Penalty weights
control.method='lbfgs';                           % Optimisation method
control.max_iter=3000;                            % Termination condition
control.parallel='ensemble';                      % Parallelisation mode

% Plots during optimisation
control.plotting={'correlation_order','coherence_order',...
                  'xy_controls','local_each_spin',...
                  'amp_controls','spectrogram'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess - smack all spins
time_axis=cumsum(control.pulse_dt);
pulse=0.05*[cos(2*pi*300*time_axis);
            ones(size(time_axis))/2];

% Run the optimisation, get normalised pulse
pulse=fminnewton(spin_system,@grape_xy,pulse);

%% Simple pulse-acquire

% Parameters
parameters.rho0=rho_init;
parameters.coil=state(spin_system,'L+','all');
parameters.pulse_op=operator(spin_system,'Ly','all');
parameters.pulse_angle=pi/4;
parameters.sweep=1000;
parameters.npoints=2048;
parameters.zerofill=4096;
parameters.axis_units='Hz';

% Simulation
fid=liquid(spin_system,@hp_acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
figure(); scale_figure([2.5 0.75]);
subplot(1,2,1); plot_1d(spin_system,real(spectrum),parameters);
xlim tight; ylim([-40 40]);

%% Optimal control pulse

% Apply power level scaling
pulse=mean(control.pwr_levels)*pulse;
pulse=mat2cell(pulse,[1 1]);

% Run a test simulation using the optimal pulse
report(spin_system,'running test simulation...');
rho=shaped_pulse_xy(spin_system,H,control.operators,pulse,...
                    control.pulse_dt,rho_init,'expv-pwc');
rho=homospoil(spin_system,rho,'destroy');
fidelity=real(rho_targ'*rho);
report(spin_system,['Re[<target|rho(T)>] = ' num2str(fidelity)]);

% New initial condition
parameters.rho0=rho;

% Simulation
fid=liquid(spin_system,@hp_acquire,parameters,'nmr');

% Apodisation
fid=apodisation(spin_system,fid,{{'exp',6}});

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
subplot(1,2,2); plot_1d(spin_system,real(spectrum),parameters);
xlim tight; ylim([-40 40]);

end

