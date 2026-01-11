% Broadband refocusing (π) pulse for liquid-state 1H NMR.
%
% Spinach implementation of the broadband refocusing example from
%  http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Task: Design a 200 µs broadband x-phase π pulse (effective Rx(pi))
% that refocuses the transverse magnetisation and inverts the longitudinal
% component:
%
%       {Sx ->  Sx,  Sy -> -Sy,  Sz -> -Sz}
%
% over an offset range of ±12.5 kHz. Robustness is enforced by an offset
% ensemble and RF power regularisation.
%
% aditya.dev@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function broadband_refocusing_pulses()

% Pulse length and discretisation
T         = 200e-6;                 % total duration, s
n_t_steps = 600;                    % number of time slices
dt        = T / n_t_steps;

% Magnetic field (Tesla)
sys.magnet   = 14.1;
sys.isotopes = {'1H'};

% Chemical shift (ppm): on-resonance reference
inter.zeeman.scalar = {0.0};

% Basis set
bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);


% Cartesian basis states
Sx = state(spin_system, 'Lx', 1);
Sy = state(spin_system, 'Ly', 1);
Sz = state(spin_system, 'Lz', 1);

% Normalise (Hilbert–Schmidt norm)
Sx = Sx / norm(full(Sx), 2);
Sy = Sy / norm(full(Sy), 2);
Sz = Sz / norm(full(Sz), 2);

% Initial and target states
rho_init = {Sx,  Sy,  Sz};
rho_targ = {Sx, -Sy, -Sz};

% RF controls and offset operator
Lx = operator(spin_system, 'Lx', 1);
Ly = operator(spin_system, 'Ly', 1);
Lz = operator(spin_system, 'Lz', 1);

% Drift Hamiltonian
H  = hamiltonian(assume(spin_system, 'nmr'));

% Control data structure
control.drifts    = {{H}};                                 % Drift Hamiltonian
control.operators = {Lx, Ly};                              % Control operators
control.off_ops   = {Lz};                                  % Offset operator
control.offsets   = {2*pi*linspace(-12.5e3, 12.5e3, 101)}; % as per the paper
control.rho_init  = rho_init;                              % Initial states
control.rho_targ  = rho_targ;                              % Target states
control.pulse_dt  = dt * ones(1, n_t_steps);               % as per the paper
control.pwr_levels = 2*pi*30e3;                            % as per the paper
control.penalties = {'NS', 'SNS'};                         % Penalties
control.p_weights = [0.01 100];                            % Penalty weights
control.method   = 'lbfgs';                                % Optimiser
control.max_iter = 200;                                    % Max iterations
control.parallel = 'ensemble';                             % Parallel mode

% Visual diagnostics
%control.plotting = {'phi_controls', ...
%    'xy_controls', ...
%    'amp_controls', ...
%    'robustness'};

n_channels = numel(control.operators);
% Random initial guess
guess = rand(n_channels, n_t_steps);

% Optimisation
spin_system = optimcon(spin_system, control);
xy_profile  = fminnewton(spin_system, @grape_xy, guess);

% Convert normalised waveform to physical rad/s controls
CLx = control.pwr_levels*xy_profile(1,:);
CLy = control.pwr_levels*xy_profile(2,:);

% Bookshelf initial/target stacks for a multi-target check
rho0 = [Sx,  Sy,  Sz];
rtg  = [Sx,  -Sy,  -Sz];

% Offset grid 
test_offsets = linspace(-12.5e3, 12.5e3, 501);

avg_eff = zeros(size(test_offsets));
eff_Sx  = zeros(size(test_offsets));
eff_Sy  = zeros(size(test_offsets));
eff_Sz  = zeros(size(test_offsets));

for k = 1:numel(test_offsets)
    w = 2*pi*test_offsets(k);           % rad/s
    Hk = H + w*Lz;                      % include offset in drift
    rho_k = shaped_pulse_xy(spin_system, Hk, control.operators, {CLx, CLy}, ...
        control.pulse_dt, rho0, 'expv-pwc');

    ov = diag(rho0' * rho_k);            % overlaps for each column
    eff_Sx(k) = real(ov(1));
    eff_Sy(k) = real(ov(2));
    eff_Sz(k) = real(ov(3));
    avg_eff(k)= mean([eff_Sx(k), eff_Sy(k), eff_Sz(k)]);
end

figure(1); clf;

subplot(2,2,1);
plot(test_offsets/1e3, eff_Sx, '-'); grid on;
xlabel('Offset (kHz)'); ylabel('Overlap');
title('Sx \rightarrow Sx');
ylim([-1.05 1.05]);

subplot(2,2,2);
plot(test_offsets/1e3, eff_Sy, '-'); grid on;
xlabel('Offset (kHz)'); ylabel('Overlap');
title('Sy \rightarrow -Sy');
ylim([-1.05 1.05]);

subplot(2,2,3);
plot(test_offsets/1e3, eff_Sz, '-'); grid on;
xlabel('Offset (kHz)'); ylabel('Overlap');
title('Sz \rightarrow -Sz');
ylim([-1.05 1.05]);

subplot(2,2,4);
plot(test_offsets/1e3, avg_eff, '--'); grid on;
xlabel('Offset (kHz)'); ylabel('Overlap');
title('Average');
ylim([-1.05 1.05]);

end
