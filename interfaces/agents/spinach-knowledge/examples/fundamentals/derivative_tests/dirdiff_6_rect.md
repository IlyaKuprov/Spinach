# examples/fundamentals/derivative_tests/dirdiff_6_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_6_rect.m`
- Signature: `dirdiff_6_rect()`
- Total lines: 85

## Purpose

GRAPE Hessian test against finite-differenced gradients.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 11-12: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 14-15: Build the derivative-test system; implemented by `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 17-18: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 34-35: Random guess and finite diff increment; implemented by `guess=randn(2,5)/3; h=sqrt(eps('double'))`.
- Lines 37-38: Call GRAPE and request analytical Hessian; implemented by `[~,~,~,hess_anl]=grape_xy(guess,spin_system)`.
- Lines 41-42: Leftmost Hessian column; implemented by `wave_forw=guess; wave_forw(1)=wave_forw(1)+h`.
- Lines 55-56: Rightmost Hessian column; implemented by `wave_forw=guess; wave_forw(end)=wave_forw(end)+h`.
- Lines 69-70: Middle Hessian column; implemented by `wave_forw=guess; wave_forw(3)=wave_forw(3)+h`.

### Control flow inferred from the code

- Line 12: `for` loop over `n=1:numel(formalisms)`.
- Line 49: conditional branch on `norm(hess_anl(:,1)-hess_num,1)<1e-6*norm(hess_num,1)`.
- Line 63: conditional branch on `norm(hess_anl(:,end)-hess_num,1)<1e-6*norm(hess_num,1)`.
- Line 77: conditional branch on `norm(hess_anl(:,3)-hess_num,1)<1e-6*norm(hess_num,1)`.

### Key state/data transformations

- Lines 9: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 15: computes `[spin_system,Sx,Sy,Sz,Lx,Ly,H]` using `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 18: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 19: computes `control.channels` using `control.channels=[1;1]`.
- Lines 20: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 21: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 22: computes `control.rho_init` using `control.rho_init={ Sx Sy Sz}`.
- Lines 23: computes `control.rho_targ` using `control.rho_targ={-Sz Sy Sx}`.
- Lines 24: computes `control.pwr_levels` using `control.pwr_levels=2*pi*linspace(50e3,70e3,10)`.
- Lines 25: computes `control.max_iter` using `control.max_iter=1000`.
- Lines 26: computes `control.plotting` using `control.plotting={}`.
- Lines 27: computes `control.integrator` using `control.integrator='rectangle'`.
- Lines 28: computes `control.pulse_dt` using `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 29: computes `control.method` using `control.method='newton'`.
- Lines 32: computes `spin_system` using `spin_system=optimcon(spin_system,control)`.
- Lines 35: computes `guess` using `guess=randn(2,5)/3; h=sqrt(eps('double'))`.
- Lines 38: computes `[~,~,~,hess_anl]` using `[~,~,~,hess_anl]=grape_xy(guess,spin_system)`.
- Lines 39: computes `hess_anl` using `hess_anl=squeeze(hess_anl(:,:,1))`.

## Implementation structure

- GRAPE Hessian test against finite-differenced gradients.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Spinach housekeeping
- Random guess and finite diff increment
- Call GRAPE and request analytical Hessian
- Leftmost Hessian column
- Rightmost Hessian column
- Middle Hessian column

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `optimcon()`, `eps()`, `grape_xy()`, `squeeze()`, `hess_anl()`, `wave_forw()`, `wave_back()`, `grad_forw()`, `grad_back()`.
