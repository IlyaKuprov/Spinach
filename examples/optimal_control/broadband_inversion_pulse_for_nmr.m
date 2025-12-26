function broadband_inversion_pulse_for_nmr()
% Broadband inversion pulse design for liquid-state NMR.
%
% This script reproduces, using Spinach, the second SIMPSON optimal control
% example from:
%
%   Z. Tošner, C. Kehlet, N. Khaneja, S.J. Glaser, N.C. Nielsen,
%   "Optimal control in NMR spectroscopy: Numerical implementation
%    in SIMPSON", J. Magn. Reson. 197 (2009) 120–134.
%
%        http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Physical problem:
% -----------------
% A single-spin (1H) system is considered in the rotating frame, on
% resonance (no isotropic chemical shift). The goal is to design a
% broadband inversion pulse that performs:
%
%       I_z  →  -I_z
%
% uniformly over a frequency offset range of ±50 kHz, under a fixed pulse
% duration T = 600 µs, discretised into 600 time steps of 1 µs.
%
% The control fields are Cartesian RF components (Lx, Ly) on 1H, and
% robustness is enforced by an offset ensemble and RF-power penalties.
%
% Authors:
%   Aditya Dev  (aditya.dev@weizmann.ac.il)
%   Ilya Kuprov (ilya.kuprov@weizmann.ac.il)

%==========================================================================
% 1. Spin system specification
%==========================================================================

% Total pulse duration (s) and time discretisation
T         = 600e-6;                  % 600 µs total duration
n_t_steps = 600;                     % 600 time steps
dt        = T / n_t_steps;           % 1 µs per step

% Static magnetic field (Tesla)
sys.magnet   = 14.1;
sys.isotopes = {'1H'};

% Isotropic chemical shift (ppm) – on resonance
inter.zeeman.scalar = {0.0};

% Basis set
bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%==========================================================================
% 2. Initial and target states
%==========================================================================

Sx = state(spin_system, 'Lx', 1);
Sy = state(spin_system, 'Ly', 1);
Sz = state(spin_system, 'Lz', 1);

Sx = Sx / norm(full(Sx), 2);
Sy = Sy / norm(full(Sy), 2);
Sz = Sz / norm(full(Sz), 2);

% State-to-state transfer: Sz -> -Sz (longitudinal inversion)
rho_init =  Sz;
rho_targ = -Sz;

%==========================================================================
% 3. Drift and control operators
%==========================================================================

LxH = operator(spin_system, 'Lx', 1);
LyH = operator(spin_system, 'Ly', 1);

% Offset operator (z-component) used to generate the offset ensemble
LzH = operator(spin_system, 'Lz', 1);

% Drift Hamiltonian (liquid-state NMR context)
H = hamiltonian(assume(spin_system, 'nmr'));

%==========================================================================
% 4. Optimal control setup
%==========================================================================

% Drift and controls
control.drifts    = {{H}};
control.operators = {LxH, LyH};

% Offset ensemble
% The SIMPSON example considers two cases:
%   (a) offsets separated by 10 kHz (11 points)
%   (b) offsets separated by  1 kHz (101 points)
% Here we consider case (a); for case (b), change 11 -> 101.
control.off_ops = {LzH};
control.offsets = {2*pi*linspace(-50e3, 50e3, 101)};   % rad/s

% Objective
control.rho_init = {rho_init};
control.rho_targ = {rho_targ};

% Time grid
control.pulse_dt = dt * ones(1, n_t_steps);

% RF power specification
% The SIMPSON example reports an RMS RF amplitude of ~13 kHz with peak
% values reaching ~30 kHz; a separate run enforces a stricter 10 kHz
% amplitude cap. Here we keep the RF scale at 13 kHz.
control.pwr_levels = 2*pi * 13e3;

% Amplitude profile (grape_xy optimises Cartesian amplitudes)
n_channels         = numel(control.operators);
control.amplitudes = ones(n_channels, n_t_steps);

% Penalties and optimiser
% 'NS' penalises accumulated RF energy.
% 'SNS' penalises overspill beyond the nominal amplitude scale; used here as
% a soft proxy for a hard RF bound (Spinach does not enforce hard bounds).
control.penalties = {'NS', 'SNS'};
control.p_weights = [0.1 100];

control.method   = 'lbfgs';
control.max_iter = 200;
control.parallel = 'ensemble';

% Plotting options (enable if needed)
% control.plotting = {'xy_controls','amp_controls','robustness'};

%==========================================================================
% 5. Initial guess
%==========================================================================

guess = randn(n_channels, n_t_steps);

%==========================================================================
% 6. Optimisation
%==========================================================================

spin_system = optimcon(spin_system, control);
xy_profile  = fminnewton(spin_system, @grape_xy, guess);

end
