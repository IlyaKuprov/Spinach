% Optimal control pulse starting with Lz and populating the
% Ly state on 87Rb in a quadrupolar rubidium system under
% magic angle spinning. A phase-modulated pulse is produced.
% 
% Calculation time: hours.
%
% i.kuprov@soton.ac.uk
% m.carravetta@soton.ac.uk

function mas_powder_control()

% System specification
sys.magnet=9.413;         % magnet
sys.isotopes={'87Rb'};    % spins

% Quadrupolar coupling
inter.coupling.matrix{1,1}=eeqq2nqi(1.68e6,0.2,3/2,[0 0 0]);

% Basis set and formalism
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% MAS experiment parameters
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];    % spinning axis
parameters.rate=-20e3;                      % spinning rate, Bruker direction
parameters.grid='rep_2ang_100pts_sph';      % powder grid
parameters.max_rank=8;                      % spinner rank
parameters.spins={'87Rb'};                  % spins involved
parameters.rframes={{'87Rb',3}};            % rotating frame order
parameters.verbose=0;                       % diagnostic output

% Drift Liouvillians and classical subspace dimension for the ensemble
[control.drifts,spc_dim]=drifts(spin_system,@singlerot,parameters,'qnmr');

% Initial state - Lz
rho_init=state(spin_system,'Lz','87Rb');
space_part=ones(spc_dim,1);
rho_init=kron(space_part,rho_init);
rho_init=rho_init/norm(rho_init,2);
control.rho_init={rho_init};

% Target state - Ly
rho_targ=state(spin_system,'Ly','87Rb');
rho_targ=kron(space_part,rho_targ);
rho_targ=rho_targ/norm(rho_targ);
control.rho_targ={rho_targ};

% Control operators
Lx=operator(spin_system,'Lx','87Rb');
Ly=operator(spin_system,'Ly','87Rb');  
space_part=speye(spc_dim,spc_dim);
Lx=kron(space_part,Lx);
Ly=kron(space_part,Ly);

% Offset operator
Lz=operator(spin_system,'Lz','87Rb'); 
Lz=kron(space_part,Lz);

% Control parameters
control.operators={Lx,Ly};                          % Control operators
control.pulse_dt=0.5e-6*ones(1,100);                % Slice up one rotor period
control.pwr_levels=2*pi*[110 120 130]*1e3/sqrt(2);  % Power per channel
control.offsets={linspace(-1e3,+1e3,5)};            % Offset distribution
control.off_ops={Lz};                               % Offset operator
control.method='lbfgs';                             % Optimisation method
control.max_iter=100;                               % Termination condition
control.amplitudes=ones(1,100);                     % Amplitude profile
control.parallel='ensemble';                        % Parallelisation

% Plotting options
control.plotting={'phi_controls','robustness',...
                  'xy_controls','spectrogram'};

% Initial guess
guess=pi/2*ones(1,100); 

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run the optimisation
pulse=fminnewton(spin_system,@grape_phase,guess);

% Construct Cartesian waveform
pulse_x=mean(control.pwr_levels)*cos(pulse);
pulse_y=mean(control.pwr_levels)*sin(pulse);

% Average over the systems
parfor n=1:numel(control.drifts)
    
    % Run the pulse
    rho=shaped_pulse_xy(spin_system,control.drifts{n}{1},control.operators,...
                        {pulse_x,pulse_y},control.pulse_dt,rho_init,'expv-pwc'); %#ok<PFBNS>
    fid(n)=real(rho_targ'*rho);
    
end

% Report fidelity
disp(['Average fidelity ' num2str(mean(fid))]);

end

