# examples/fundamentals/derivative_tests/dirdiff_5_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_5_rect.m`
- Signature: `dirdiff_5_rect()`
- Total lines: 71

## Purpose

GRAPE Hessian internal consistency test: Newton against Goodwin algorithm.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 12-13: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 15-16: Build the derivative-test system; implemented by `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 18-19: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 31-32: Pick initial guess, phase-modulated GRAPE; implemented by `control.amplitudes=ones(1,5); guess=randn(1,5)/3`.
- Lines 34-35: Get Newton Hessian; implemented by `control.method='newton'`.
- Lines 39-40: Get Goodwin Hessian; implemented by `control.method='goodwin'`.
- Lines 44-45: Pick initial guess, XY-modulated GRAPE; implemented by `guess=randn(2,5)/3`.
- Lines 57-58: Run the comparisons; implemented by `if norm(newton_hess_ph(:)-goodwin_hess_ph(:),1)>1e-6*norm(newton_hess_ph(:),1)`.

### Control flow inferred from the code

- Line 13: `for` loop over `n=1:numel(formalisms)`.
- Line 58: conditional branch on `norm(newton_hess_ph(:)-goodwin_hess_ph(:),1)>1e-6*norm(newton_hess_ph(:),1)`.
- Line 63: conditional branch on `norm(newton_hess_xy(:)-goodwin_hess_xy(:),1)>1e-6*norm(newton_hess_xy(:),1)`.

### Key state/data transformations

- Lines 10: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 16: computes `[spin_system,Sx,Sy,Sz,Lx,Ly,H]` using `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 19: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 20: computes `control.channels` using `control.channels=[1;1]`.
- Lines 21: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 22: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 23: computes `control.rho_init` using `control.rho_init={ Sx Sy Sz}`.
- Lines 24: computes `control.rho_targ` using `control.rho_targ={-Sz Sy Sx}`.
- Lines 25: computes `control.pwr_levels` using `control.pwr_levels=2*pi*linspace(50e3,70e3,10)`.
- Lines 26: computes `control.max_iter` using `control.max_iter=1000`.
- Lines 27: computes `control.plotting` using `control.plotting={}`.
- Lines 28: computes `control.integrator` using `control.integrator='rectangle'`.
- Lines 29: computes `control.pulse_dt` using `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 32: computes `control.amplitudes` using `control.amplitudes=ones(1,5); guess=randn(1,5)/3`.
- Lines 35: computes `control.method` using `control.method='newton'`.
- Lines 36: computes `spin_system` using `spin_system=optimcon(spin_system,control)`.
- Lines 37: computes `[~,~,~,newton_hess_ph]` using `[~,~,~,newton_hess_ph]=grape_phase(guess,spin_system)`.
- Lines 42: computes `[~,~,~,goodwin_hess_ph]` using `[~,~,~,goodwin_hess_ph]=grape_phase(guess,spin_system)`.

## Implementation structure

- GRAPE Hessian internal consistency test: Newton
- against Goodwin algorithm.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Pick initial guess, phase-modulated GRAPE
- Get Newton Hessian
- Get Goodwin Hessian
- Pick initial guess, XY-modulated GRAPE
- Run the comparisons

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `optimcon()`, `grape_phase()`, `grape_xy()`, `newton_hess_ph()`, `goodwin_hess_ph()`, `newton_hess_xy()`, `goodwin_hess_xy()`.
