% Transmitter offset selective excitation described in Glaser group 
% paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified
% transmitter offset intervals have magnetisation arriving into us-
% er specified states. The pulse is phase-modulated.
%
% Calculation time: minutes.
%
% i.kuprov@soton.ac.uk

function pattern_pulse_2()

% Magnetic field
sys.magnet=28.18;

% Single carbon spin
sys.isotopes{1}='13C';

% Put the spin at 0 ppm
inter.zeeman.scalar{1}=0;

% No approximations
bas.formalism='sphten-liouv';
bas.approximation='none';

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get pertinent spin states
Sx=state(spin_system,'Lx','13C'); Sx=Sx/norm(Sx,2);
Sz=state(spin_system,'Lz','13C'); Sz=Sz/norm(Sz,2);

% Get pertinent control operators
Lx=operator(spin_system,'Lx','13C');
Ly=operator(spin_system,'Ly','13C');
Lz=operator(spin_system,'Lz','13C');

% Get the drift Hamiltonian (zero here)
H=hamiltonian(assume(spin_system,'nmr'));

% Transmitter offset range in Hz
toff_range=linspace(-4e3,4e3,128);

% Excitation pattern
X=ones(1,128);  X(1:20)=0; X(109:128)=0; X(54:73)=0;
Z=zeros(1,128); Z(1:20)=1; Z(109:128)=1; Z(54:73)=1;
figure(2); plot(toff_range,[X; Z]); kgrid;
kxlabel('Transmitter offset, Hz'); 
klegend({'X (target)','Z (target)'});
xlim tight; ylim([-0.1 1.1]); drawnow; figure(1); 

% Patterened target state array
rho_targ=X.*repmat(Sx,1,128)+Z.*repmat(Sz,1,128);

% Define control parameters
control.drifts={{H}};                     % Drift
control.operators={Lx,Ly};                % Controls
control.off_ops={Lz};                     % Offset operator
control.offsets={toff_range};             % Transmitter offsets
control.rho_init=repmat({Sz},1,128);      % Starting states
control.rho_targ=num2cell(rho_targ,1);    % Target states
control.pulse_dt=2e-5*ones(1,300);        % Pulse interval grid
control.pwr_levels=2*pi*2000;             % Power level
control.amplitudes=ones(1,300);           % Amplitude profile
control.method='lbfgs';                   % Optimisation method
control.max_iter=200;                     % Termination condition
control.parallel='ensemble';              % Parallelisation
control.ens_corrs={'rho_ens'};            % Own state pair for each B1

% Plotting options
control.plotting={'phi_controls','xy_controls',...
                  'robustness','spectrogram'};

% Initial guess
guess=(pi/4)*ones(1,300);

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run LBFGS GRAPE pulse optimisation
phi_profile=fminnewton(spin_system,@grape_phase,guess);

% Loop over the offsets
X=zeros(1,128); Z=zeros(1,128);
parfor n=1:numel(toff_range) %#ok<*PFBNS>

    % Get Cartesian components of the pulse
    amp_profile=control.pwr_levels*control.amplitudes; 
    [CLx,CLy]=polar2cartesian(amp_profile,phi_profile);

    % Simulate the optimised pulse
    rho_init=control.rho_init{n};
    rho=shaped_pulse_xy(spin_system,H+2*pi*toff_range(n)*control.off_ops{1},...
                        {Lx,Ly},{CLx,CLy},control.pulse_dt,rho_init,'expv-pwc');

    % Project out the destination states
    X(n)=real(Sx'*rho); Z(n)=real(Sz'*rho);

end

% Update the plot
figure(2); hold on; plot(toff_range,[X; Z],'o');
klegend({'X (target)','Z (target)','X (result)','Z (result)'});

end