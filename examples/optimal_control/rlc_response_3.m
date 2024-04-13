% Probe circuit response effect on the accuracy of the deu-
% terium pre-phasing pulse designed to set deuterium magne-
% tisation in a -CD3 group of alanine up for rephasing 100
% microseconds after the pulse is finished.
%
% The system is assumed to be a powder (100 orientations) 
% with a B1 distribution (from 40 to 60 kHz per channel).
%
% Calculation time: minutes
%
% i.kuprov@soton.ac.uk
% marina.carravetta@soton.ac.uk
% uluk.rasulov@soton.ac.uk

function rlc_response_3()

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
parameters.grid='rep_2ang_200pts_sph';
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

% Get control operators
Lx=operator(spin_system,'Lx','2H');
Ly=operator(spin_system,'Ly','2H');

% Define control parameters
control.operators={Lx,Ly};                     % Controls
control.rho_init={rho_init};                   % Starting state
control.rho_targ={rho_targ};                   % Destination state
control.pwr_levels=2*pi*35e3;                  % Power per channel
control.pulse_dt=6.5e-6*ones(1,24);            % Slice duration grid
control.method='lbfgs';                        % Optimisation method
control.max_iter=1000;                         % Maximum iterations
control.dead_time=100e-6;                      % Dead time
control.penalties={'SNS'};                     % Penalty types
control.p_weights=100;                         % Penalty weights
control.integrator='trapezium';                % Piecewise-constant
control.parallel='ensemble';                   % Parallelisation

% Freeze the edges
control.freeze=false(2,25);                   
control.freeze(:,1:2)=true;
control.freeze(:,(end-1):end)=true;

% Plotting options
control.plotting={'xy_controls','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Waveform guess with zeros at the edges
guess=randn(2,25)/10; guess(:,1:2)=0; 
guess(:,(end-1):end)=0;

% Run the optimisation
pulse_profile=fminnewton(spin_system,@grape_xy,guess);

% Get Cartesian components
CLx=mean(control.pwr_levels)*pulse_profile(1,:);
CLy=mean(control.pwr_levels)*pulse_profile(2,:);

% Apply RLC response distortion by a probe circuit 
% with Q = 200 and time shift compensation
figure(); restrans(CLx',CLy',control.pulse_dt(1),...
                   sys.magnet*spin('2H'),200,'pwl_tsc',100);

end

