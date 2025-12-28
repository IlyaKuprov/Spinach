function broadband_refocusing_pulses()
% Broadband refocusing (π) pulse for liquid-state 1H NMR.
%
% Spinach implementation of the broadband refocusing example from:
%
%   Z. Tošner, C. Kehlet, N. Khaneja, S.J. Glaser, N.C. Nielsen, J. Magn.
%   Reson. 197 (2009) 120–134.
%
%             http://dx.doi.org/10.1016/j.jmr.2008.11.020
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
% Authors:
%   Aditya Dev  (aditya.dev@weizmann.ac.il) 
%   Ilya Kuprov (ilya.kuprov@weizmann.ac.il)

%---------------------------------------------------------------------------%
% Spin system
%---------------------------------------------------------------------------%

% Pulse length and discretisation
T         = 200e-6;                 % total duration, s
n_t_steps = 600;                    % number of time slices
dt        = T / n_t_steps;          % slice duration

% Static field and isotope
sys.magnet   = 14.1;                % Tesla
sys.isotopes = {'1H'};

% Chemical shift (ppm): zero => on-resonance reference
inter.zeeman.scalar = {0.0};

% Basis set: full Liouville space in spherical tensor formalism
bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

% Spinach housekeeping
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%---------------------------------------------------------------------------%
% States and operators
%---------------------------------------------------------------------------%

% Cartesian basis states for verification / multi-target optimisation
Sx = state(spin_system, 'Lx', 1);
Sy = state(spin_system, 'Ly', 1);
Sz = state(spin_system, 'Lz', 1);

% Normalise (Hilbert–Schmidt norm)
Sx = Sx / norm(full(Sx), 2);
Sy = Sy / norm(full(Sy), 2);
Sz = Sz / norm(full(Sz), 2);

% Multi-target refocusing objective: Rx(pi) leaves Sx invariant and flips
% the signs of Sy and Sz
rho_init = {Sx,  Sy,  Sz};
rho_targ = {Sx, -Sy, -Sz};

% RF controls (Cartesian components) and offset operator
Lx = operator(spin_system, 'Lx', 1);
Ly = operator(spin_system, 'Ly', 1);
Lz = operator(spin_system, 'Lz', 1);

% Drift Hamiltonian (liquid-state NMR context)
H  = hamiltonian(assume(spin_system, 'nmr'));

%---------------------------------------------------------------------------%
% Optimal control settings
%---------------------------------------------------------------------------%

% Drift and control operators
control.drifts    = {{H}};
control.operators = {Lx, Ly};

% Offset ensemble: robustness over ±12.5 kHz (25 kHz bandwidth total)
control.off_ops   = {Lz};
control.offsets   = {2*pi*linspace(-12.5e3, 12.5e3, 100)};   % rad/s

% Multi-target specification
control.rho_init  = rho_init;
control.rho_targ  = rho_targ;

% Piecewise-constant time grid
control.pulse_dt  = dt * ones(1, n_t_steps);

% RF scaling (nominal amplitude scale, ~15 kHz)
control.pwr_levels = 2*pi*15e3;

% Regularisation:
%   NS  - suppresses accumulated RF energy 
%   SNS - penalises control spillout beyond nominal amplitude bounds
control.penalties = {'NS', 'SNS'};
control.p_weights = [0.01 100];

% Optimiser configuration
control.method   = 'lbfgs';
control.max_iter = 200;
control.parallel = 'ensemble';

% Visual diagnostics during optimisation
control.plotting = {'phi_controls', ...
    'xy_controls', ...
    'amp_controls', ...
    'robustness'};

%---------------------------------------------------------------------------%
% Initial guess and optimisation
%---------------------------------------------------------------------------%

% Two controls (Lx, Ly) => two waveform rows
n_channels = numel(control.operators);

% Initial guess waveform (Cartesian components)
guess = rand(n_channels, n_t_steps);

% Flat initial amplitude mask
control.amplitudes = ones(n_channels, n_t_steps);

% Spinach housekeeping and run optimisation
spin_system = optimcon(spin_system, control);
xy_profile  = fminnewton(spin_system, @grape_xy, guess); 

end
