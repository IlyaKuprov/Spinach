# examples/fundamentals/derivative_tests/dirdiff_3_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_3_rect.m`
- Signature: `dirdiff_3_rect()`
- Total lines: 82

## Purpose

Directional derivative test for the Cartesian GRAPE module, rectangles integrator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 12-13: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 15-16: Build the derivative-test system; implemented by `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 18-19: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 31-32: Set the interval grid; implemented by `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 37-38: Random guess and finite diff increment; implemented by `guess=randn(2,5)/3; h=sqrt(eps('double'))`.
- Lines 40-41: Call GRAPE and request analytical gradient; implemented by `[~,~,grad_anl]=grape_xy(guess,spin_system)`.
- Lines 44-45: Left waveform edge; implemented by `wave_forw=guess; wave_forw(1)=wave_forw(1)+h`.
- Lines 56-57: Right waveform edge; implemented by `wave_forw=guess; wave_forw(end)=wave_forw(end)+h`.
- Lines 68-69: Waveform midpoint; implemented by `wave_forw=guess; wave_forw(3)=wave_forw(3)+h`.

### Control flow inferred from the code

- Line 13: `for` loop over `n=1:numel(formalisms)`.
- Line 50: conditional branch on `abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6`.
- Line 62: conditional branch on `abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6`.
- Line 74: conditional branch on `abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6`.

### Key state/data transformations

- Lines 10: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 16: computes `[spin_system,Sx,Sy,Sz,Lx,Ly,H]` using `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 19: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 20: computes `control.channels` using `control.channels=[1;1]`.
- Lines 21: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 22: computes `control.operators` using `control.operators={Lx Ly}`.
- Lines 23: computes `control.rho_init` using `control.rho_init={ Sx Sy Sz}`.
- Lines 24: computes `control.rho_targ` using `control.rho_targ={-Sz Sy Sx}`.
- Lines 25: computes `control.pwr_levels` using `control.pwr_levels=2*pi*linspace(50e3,70e3,10)`.
- Lines 26: computes `control.method` using `control.method='lbfgs'`.
- Lines 27: computes `control.max_iter` using `control.max_iter=1000`.
- Lines 28: computes `control.plotting` using `control.plotting={}`.
- Lines 29: computes `control.integrator` using `control.integrator='rectangle'`.
- Lines 32: computes `control.pulse_dt` using `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 35: computes `spin_system` using `spin_system=optimcon(spin_system,control)`.
- Lines 38: computes `guess` using `guess=randn(2,5)/3; h=sqrt(eps('double'))`.
- Lines 41: computes `[~,~,grad_anl]` using `[~,~,grad_anl]=grape_xy(guess,spin_system)`.
- Lines 42: computes `grad_anl` using `grad_anl=squeeze(grad_anl(:,:,1))`.

## Implementation structure

- Directional derivative test for the Cartesian GRAPE
- module, rectangles integrator.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Set the interval grid
- Spinach housekeeping
- Random guess and finite diff increment
- Call GRAPE and request analytical gradient
- Left waveform edge
- Right waveform edge

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `optimcon()`, `eps()`, `grape_xy()`, `squeeze()`, `grad_anl()`, `wave_forw()`, `wave_back()`, `fid_forw()`, `fid_back()`.
