% Coherence transfer through scalar J-coupling in a two-spin system.
%
% This script implements, in Spinach, the first SIMPSON optimal control
% example from http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% A heteronuclear two-spin system (1H–13C) with an isotropic scalar
% coupling J, with both nuclei set on resonance (no chemical shift
% offsets). The goal is to transfer transverse magnetisation from proton to
% carbon:
%
%       I_x(1H)  →  I_x(13C)
%
% over a fixed evolution period T = 1/J.
%
%  aditya.dev@weizmann.ac.il
%  ilya.kuprov@weizmann.ac.il

function coherence_transfer_through_J_coup()

J = 140;  % Scalar coupling Hz
T = 1 / J; % pulse duration
n_t_steps = 150;  % Number of time steps
dt        = T / n_t_steps;

% Magnetic field (Tesla)
sys.magnet = 14.1;
sys.isotopes = {'1H','13C'};

% Isotropic chemical shifts
inter.zeeman.scalar = {0.0, 0.0};

% Isotropic scalar coupling
inter.coupling.scalar        = cell(2,2);
inter.coupling.scalar{1,2}   = J;     % Hz

% Basis set
bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

% Initial state: Lx on proton (spin 1)
rho_init = state(spin_system, {'Lx'}, {1});
rho_init = rho_init / norm(full(rho_init), 2);

% Target state: Lx on carbon (spin 2)
rho_targ = state(spin_system, {'Lx'}, {2});
rho_targ = rho_targ / norm(full(rho_targ), 2);

% RF control operators on 1H
LxH = operator(spin_system, 'Lx', 1);
LyH = operator(spin_system, 'Ly', 1);

% RF control operators on 13C
LxC = operator(spin_system, 'Lx', 2);
LyC = operator(spin_system, 'Ly', 2);

% Drift Hamiltonian
H = hamiltonian(assume(spin_system, 'nmr'));


% Control data structure
control.drifts = {{H}};                           % Drift Hamiltonian
control.operators = {LxH, LyH, LxC, LyC};         % Control operators
control.rho_init = {rho_init};                    % Initial state
control.rho_targ = {rho_targ};                    % Target state
control.pulse_dt = dt * ones(1, n_t_steps);       % Time discretisation
control.method   = 'lbfgs';                       % Optimiser
control.max_iter = 200;                           % Max iterations
control.fidelity = 'real';
control.parallel = 'ensemble';

% Case (a) is without any penalty check Case (b) imposes the penalty on
% accumulated RF power We implement case (b) here.
control.penalties = {'NS'};                          % Penalties
control.p_weights = 0.01;                            % Penalty weight

% RF power levels In actual implementation, the SIMPSON example takes 30
% random rf amplitudes from 0 to 1000 or 1500 Hz range for H and C channel
% respectively Here we choose 30 random values between 0 to 1000 Hz over
% all channels.
control.pwr_levels = 2*pi * 1000 * rand(1,30);

n_channels = numel(control.operators);
guess      = rand(n_channels, n_t_steps); % Random initial guess

% Optimisation
spin_system = optimcon(spin_system, control);
xy_profile  = fminnewton(spin_system, @grape_xy, guess);

t = ((1:n_t_steps)-0.5)*dt;  % time axis (s)

% mean rf
rf_nom = mean(control.pwr_levels);    % rad/s

% Convert normalised waveform
uLxH = rf_nom*xy_profile(1,:);
uLyH = rf_nom*xy_profile(2,:);
uLxC = rf_nom*xy_profile(3,:);
uLyC = rf_nom*xy_profile(4,:);

figure(1); clf;

subplot(2,1,1);
plot(t*1e3, uLxH/(2*pi*1e3), t*1e3, uLyH/(2*pi*1e3), 'LineWidth', 1); grid on;
xlabel('time (ms)'); ylabel('^{1}H RF (kHz)');
title('Optimised pulse (Cartesian components)');
legend({'Lx(^{1}H)','Ly(^{1}H)'}, 'Location','best');

subplot(2,1,2);
plot(t*1e3, uLxC/(2*pi*1e3), t*1e3, uLyC/(2*pi*1e3), 'LineWidth', 1); grid on;
xlabel('time (ms)'); ylabel('^{13}C RF (kHz)');
legend({'Lx(^{13}C)','Ly(^{13}C)'}, 'Location','best');

% Robustness check over the RF ensemble points

pwr_dense = linspace(min(control.pwr_levels), max(control.pwr_levels), 101).';   % rad/s
eff_dense = zeros(numel(pwr_dense),1);

for k = 1:numel(pwr_dense)

    u = { pwr_dense(k)*xy_profile(1,:), pwr_dense(k)*xy_profile(2,:), ...
        pwr_dense(k)*xy_profile(3,:), pwr_dense(k)*xy_profile(4,:) };

    rho_f = shaped_pulse_xy(spin_system, H, control.operators, u, ...
        control.pulse_dt, rho_init, 'expv-pwc');

    eff_dense(k) = real(full(rho_targ' * rho_f));
end

figure(2); clf;
plot(pwr_dense/(2*pi), eff_dense, '-'); grid on;
xlabel('RF scaling (Hz)'); ylabel('Re(<rho_{targ} | rho(T)>)');
title('Transfer efficiency vs RF scaling (dense scan)');
ylim([-0.05 1.05]);
end
