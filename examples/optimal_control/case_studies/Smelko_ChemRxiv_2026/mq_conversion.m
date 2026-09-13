% Optimal control design of the multiple-quantum conversion pulse
% of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi-
% nach, the conversion pulse optimisation from
%
%           https://doi.org/10.26434/chemrxiv.15008427
%
% A single 27Al nucleus with the quadrupolar coupling and the shi-
% elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10
% ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz
% magnet. The quadrupolar interaction is taken to second order in
% the rotating frame, and the powder average runs over 100 crystal-
% lite orientations at 20 initial rotor phases each. The pulse is
% three rotor periods (240 us) long in 0.5 us slices, the controls
% are Cartesian, and the 100 kHz amplitude ceiling is enforced by
% a spillout penalty followed by clipping. The initial state is the
% Hermitian combination of the +MQ and -MQ coherences between the
% m=+3/2 and m=-3/2 levels (3Q) or between the m=+5/2 and m=-5/2
% levels (5Q), and the target is the population difference across
% the central transition, as in the paper. The resulting waveform
% is saved for the MQMAS efficiency calculation.
%
% Calculation time: hours on 128 cores.
%
% ilya.kuprov@weizmann.ac.il

function mq_conversion()

% Coherence order to convert, 3 or 5
mq_order=5;

% 400 MHz magnet
sys.magnet=2*pi*400e6/spin('1H');
sys.isotopes={'27Al'};

% Quadrupolar coupling and shielding anisotropy
inter.coupling.matrix{1,1}=eeqq2nqi(3.0e6,1.0,5/2,[0 0 0]);
inter.zeeman.eigs={[-5 -5 10]};
inter.zeeman.euler={[0 0 0]};

% Hilbert space formalism
bas.formalism='zeeman-hilb';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);
spin_system=assume(spin_system,'labframe');

% Rotor phase resolved drift Hamiltonians
parameters.spins={'27Al'};
parameters.axis=[sqrt(2/3) 0 sqrt(1/3)];
parameters.grid='rep_2ang_100pts_sph';
parameters.n_ticks=160;
parameters.n_phases=20;
parameters.n_slices=480;
control.drifts=mqmas_drifts(spin_system,parameters);

% Initial state, symmetric MQ coherence between m=+mq_order/2 and m=-mq_order/2
rho_init=zeros(6); rho_init(3.5-mq_order/2,3.5+mq_order/2)=1;
rho_init=rho_init+rho_init'; rho_init=rho_init/norm(rho_init,'fro');

% Target state, population difference across the central transition
rho_targ=diag([0 0 1 -1 0 0]); rho_targ=rho_targ/norm(rho_targ,'fro');

% Control operators
Lx=operator(spin_system,'Lx','27Al');
Ly=operator(spin_system,'Ly','27Al');

% Control parameters
control.isotopes={'27Al'};                 % Isotopes
control.channels=[1; 1];                   % Channel map
control.operators={Lx,Ly};                 % Control operators
control.rho_init={rho_init};               % Initial state
control.rho_targ={rho_targ};               % Target state
control.pulse_dt=0.5e-6*ones(1,480);       % Slice durations
control.pwr_levels=2*pi*100e3;             % Amplitude ceiling
control.penalties={'SNSA'};                % Amplitude spillout penalty
control.p_weights=100;                     % Penalty weight
control.method='lbfgs';                    % Optimisation method
control.max_iter=500;                      % Termination condition

% Plotting options
control.plotting={'amp_controls','phi_controls','spectrogram'};

% Random initial guess, amplitudes up to 10% of the ceiling
amp=0.1*rand(1,480); phi=2*pi*rand(1,480); amp(randi(480))=1;
guess=[amp.*cos(phi); amp.*sin(phi)];

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Run the optimisation
pulse=fmaxnewton(spin_system,@grape_xy,guess);

% Clip the amplitude to the ceiling
[amp,phi]=cartesian2polar(pulse(1,:),pulse(2,:)); amp=min(amp,1);
[pulse(1,:),pulse(2,:)]=polar2cartesian(amp,phi);

% Report the fidelity of the clipped pulse
[~,fidelity]=grape_xy(pulse,spin_system);
disp(['Fidelity after amplitude clipping: ' num2str(fidelity(1))]);

% Save the waveform in rad/s
pulse=control.pwr_levels*pulse; pulse_dt=control.pulse_dt;
save(['mq_conv_' num2str(mq_order) 'q.mat'],'pulse','pulse_dt');

end

