function coherence_transfer_through_J_coup()
% Coherence transfer through scalar J-coupling in a two-spin system.
%
% This script implements, in Spinach, the first SIMPSON optimal control
% example from:
%
%   Z. Tošner, C. Kehlet, N. Khaneja, S.J. Glaser, N.C. Nielsen,
%   "Optimal control in NMR spectroscopy: Numerical implementation
%    in SIMPSON", J. Magn. Reson. 197 (2009) 120–134.
%
%           http://dx.doi.org/10.1016/j.jmr.2008.11.020
%
% Physical problem:
% -----------------
% A heteronuclear two-spin system (1H–13C) with an isotropic scalar
% coupling J, with both nuclei set on resonance (no chemical shift
% offsets). The goal is to transfer transverse magnetisation from
% proton to carbon:
%
%       I_x(1H)  →  I_x(13C)
%
% over a fixed evolution period T = 1/J.
%
% Authors:
%   Aditya Dev  (aditya.dev@weizmann.ac.il)
%   Ilya Kuprov (ilya.kuprov@weizmann.ac.il)

%==========================================================================
% 1. Spin system specification
%==========================================================================

% Scalar coupling strength
J = 140;                              % Hz

% Total evolution time
T = 1 / J;                            % s

% Time discretisation
n_t_steps = 150;
dt        = T / n_t_steps;

% Static magnetic field
sys.magnet = 14.1;                    % Tesla

% Two-spin heteronuclear system
sys.isotopes = {'1H','13C'};

% Isotropic chemical shifts (on resonance)
inter.zeeman.scalar = {0.0, 0.0};     % ppm

% Isotropic scalar coupling
inter.coupling.scalar        = cell(2,2);
inter.coupling.scalar{1,2}   = J;     % Hz

% Basis set
bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%==========================================================================
% 2. Initial and target states
%==========================================================================

% Initial state: Lx on proton (spin 1)
rho_init = state(spin_system, {'Lx'}, {1});
rho_init = rho_init / norm(full(rho_init), 2);

% Target state: Lx on carbon (spin 2)
rho_targ = state(spin_system, {'Lx'}, {2});
rho_targ = rho_targ / norm(full(rho_targ), 2);

%==========================================================================
% 3. Drift and control operators
%==========================================================================

% RF control operators on 1H
LxH = operator(spin_system, 'Lx', 1);
LyH = operator(spin_system, 'Ly', 1);

% RF control operators on 13C
LxC = operator(spin_system, 'Lx', 2);
LyC = operator(spin_system, 'Ly', 2);

% Drift Hamiltonian (liquid-state NMR context)
H = hamiltonian(assume(spin_system, 'nmr'));

%==========================================================================
% 4. Optimal control setup
%==========================================================================

% Single drift Hamiltonian
control.drifts = {{H}};

% Cartesian RF controls: x/y on both spins
control.operators = {LxH, LyH, LxC, LyC};

% control.off_ops  = {LzH, LzC, JOp};
% control.offsets = {[-10 0 10], [-10 0 10], [-10 0 10]};

% State-to-state objective
control.rho_init = {rho_init};
control.rho_targ = {rho_targ};

% Time grid
control.pulse_dt = dt * ones(1, n_t_steps);

% RF power levels
% In actual implementation, the SIMPSON example
% takes 30 random rf amplftudes from
% 0 to 1000 or 1500 Hz range for H and C channel respectively
% Here we choose 30 random values between 0 to 1000 Hz over all channels.
control.pwr_levels = 2*pi * 1000 * rand(1,30);

% Optimisation method
control.method   = 'lbfgs';
control.max_iter = 200;

% RF power penalty; to keep a check on accumulated power in the system.
% Case (a) is without any penalty check, in which case we can simply
% comment out this part. Case (b) imposes the penalty on accumulated
% powers. we implement case (b) here.
control.penalties = {'NS'};
control.p_weights = 0.01;

control.fidelity = 'real';

% Parallelisation over time
control.parallel = 'time';

% Plotting options
%control.plotting = {'local_each_spin', ...
%    'phi_controls', ...
%    'xy_controls'};

%==========================================================================
% 5. Initial guess
%==========================================================================

n_channels = numel(control.operators);
guess      = rand(n_channels, n_t_steps);

%==========================================================================
% 6. Optimisation
%==========================================================================

spin_system = optimcon(spin_system, control);
xy_profile  = fminnewton(spin_system, @grape_xy, guess); 

end
