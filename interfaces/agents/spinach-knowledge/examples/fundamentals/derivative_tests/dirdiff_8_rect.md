# examples/fundamentals/derivative_tests/dirdiff_8_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_8_rect.m`
- Signature: `dirdiff_8_rect()`
- Total lines: 88

## Purpose

GRAPE phase Hessian test against finite-differenced gradients, rectangles integrator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 12-13: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 15-16: Build the derivative-test system; implemented by `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 18-19: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 31-32: Set the interval grid; implemented by `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 38-39: Random phases and finite diff increment; implemented by `guess=randn(2,5)/3; h=1e-5`.
- Lines 41-42: Call GRAPE and request analytical Hessian; implemented by `[~,~,~,hess_anl]=grape_phase(guess,spin_system)`.
- Lines 45-46: Leftmost Hessian column; implemented by `wave_forw=guess; wave_forw(1)=wave_forw(1)+h`.
- Lines 59-60: Rightmost Hessian column; implemented by `wave_forw=guess; wave_forw(end)=wave_forw(end)+h`.
- Lines 73-74: Middle Hessian column; implemented by `wave_forw=guess; wave_forw(5)=wave_forw(5)+h`.

### Control flow inferred from the code

- Line 13: `for` loop over `n=1:numel(formalisms)`.
- Line 53: conditional branch on `norm(hess_anl(:,1)-hess_num,1)<1e-5*norm(hess_num,1)`.
- Line 67: conditional branch on `norm(hess_anl(:,end)-hess_num,1)<1e-5*norm(hess_num,1)`.
- Line 81: conditional branch on `norm(hess_anl(:,5)-hess_num,1)<1e-5*norm(hess_num,1)`.

### Key state/data transformations

- Lines 10: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 16: computes `[spin_system,Sx,Sy,Sz,Lx,Ly,H]` using `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 19: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 20: computes `control.channels` using `control.channels=[1;1;1;1]`.
- Lines 21: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 22: computes `control.operators` using `control.operators={Lx,Ly,0.4*Lx+0.2*Ly,0.7*Ly-0.1*Lx}`.
- Lines 23: computes `control.rho_init` using `control.rho_init={ Sx Sy Sz}`.
- Lines 24: computes `control.rho_targ` using `control.rho_targ={-Sz Sy Sx}`.
- Lines 25: computes `control.pwr_levels` using `control.pwr_levels=2*pi*linspace(50e3,70e3,10)`.
- Lines 26: computes `control.method` using `control.method='newton'`.
- Lines 27: computes `control.max_iter` using `control.max_iter=1000`.
- Lines 28: computes `control.plotting` using `control.plotting={}`.
- Lines 29: computes `control.integrator` using `control.integrator='rectangle'`.
- Lines 32: computes `control.pulse_dt` using `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 33: computes `control.amplitudes` using `control.amplitudes=[ones(1,5); 0.8+0.1*(1:5)]`.
- Lines 36: computes `spin_system` using `spin_system=optimcon(spin_system,control)`.
- Lines 39: computes `guess` using `guess=randn(2,5)/3; h=1e-5`.
- Lines 42: computes `[~,~,~,hess_anl]` using `[~,~,~,hess_anl]=grape_phase(guess,spin_system)`.

## Implementation structure

- GRAPE phase Hessian test against finite-differenced gradients,
- rectangles integrator.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Set the interval grid
- Spinach housekeeping
- Random phases and finite diff increment
- Call GRAPE and request analytical Hessian
- Leftmost Hessian column
- Rightmost Hessian column

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `optimcon()`, `grape_phase()`, `squeeze()`, `hess_anl()`, `wave_forw()`, `wave_back()`, `grad_forw()`, `grad_back()`.
