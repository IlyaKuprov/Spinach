% Optimal control design of the central transition selective pulse
% of the z-filtered 27Al MQMAS experiment. Reproduces, using Spi-
% nach, the soft pulse optimisation from
%
%           https://doi.org/10.26434/chemrxiv.15008427
%
% A single 27Al nucleus with the quadrupolar coupling and the shi-
% elding anisotropy assumed in the paper (CQ=3.0 MHz, eta=1.0, 10
% ppm axial shielding anisotropy) is spun at 12.5 kHz in a 400 MHz
% magnet. The quadrupolar interaction is taken to second order in
% the rotating frame, and the powder average runs over 200 crystal-
% lite orientations at 80 initial rotor phases each. The pulse is
% 50 us long in 0.5 us slices, the controls are Cartesian, and the
% 10 kHz amplitude ceiling is enforced by a spillout penalty follo-
% wed by clipping. The initial state is the population difference
% across the central transition, and the target is the single-qu-
% antum coherence of the central transition, as in the paper. The
% resulting waveform is saved for the MQMAS efficiency calculation;
% the waveform supplied in this folder reached a fidelity of 0.70
% after 500 iterations, against the maximum of 1/sqrt(2) for this
% initial and target state pair.
%
% Calculation time: about an hour on 128 cores.
%
% ilya.kuprov@weizmann.ac.il

function ct_selective()

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
parameters.grid='rep_2ang_200pts_sph';
parameters.n_ticks=160;
parameters.n_phases=80;
parameters.n_slices=100;
control.drifts=mqmas_drifts(spin_system,parameters);

% Initial state, population difference across the central transition
rho_init=diag([0 0 1 -1 0 0]); rho_init=rho_init/norm(rho_init,'fro');

% Target state, single-quantum coherence of the central transition
rho_targ=zeros(6); rho_targ(3,4)=1;

% Control operators
Lx=operator(spin_system,'Lx','27Al');
Ly=operator(spin_system,'Ly','27Al');

% Control parameters
control.isotopes={'27Al'};                 % Isotopes
control.channels=[1; 1];                   % Channel map
control.operators={Lx,Ly};                 % Control operators
control.rho_init={rho_init};               % Initial state
control.rho_targ={rho_targ};               % Target state
control.pulse_dt=0.5e-6*ones(1,100);       % Slice durations
control.pwr_levels=2*pi*10e3;              % Amplitude ceiling
control.penalties={'SNSA'};                % Amplitude spillout penalty
control.p_weights=100;                     % Penalty weight
control.method='lbfgs';                    % Optimisation method
control.max_iter=500;                      % Termination condition

% Plotting options
control.plotting={'amp_controls','phi_controls','spectrogram'};

% Random initial guess, amplitudes up to 10% of the ceiling
amp=0.1*rand(1,100); phi=2*pi*rand(1,100); amp(randi(100))=1;
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
save('ct_pulse.mat','pulse','pulse_dt');

end

