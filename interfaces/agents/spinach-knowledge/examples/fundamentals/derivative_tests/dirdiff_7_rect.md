# examples/fundamentals/derivative_tests/dirdiff_7_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_7_rect.m`
- Signature: `dirdiff_7_rect()`
- Total lines: 88

## Purpose

Directional derivative test for the phase-modulated GRAPE module, with an ensemble of chained linear and non-linear instrumental distortions of the waveform.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Formalisms to test; implemented by `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 13-14: Loop over formalisms; implemented by `for n=1:numel(formalisms)`.
- Lines 16-17: Build the derivative-test system; implemented by `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 19-20: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 32-33: Define an ensemble of distortion chains; implemented by `control.distortion={@(w)firf(w,[0.9 0.1i]), @(w)spf(w,0.2), @(w)szf(w,0.2), @(w)amp_root(w,2*pi*20e3,4)`.
- Lines 36-37: Set the interval grid; implemented by `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 43-44: Random phases and finite diff increment; implemented by `guess=randn(1,5)/3; h=sqrt(eps('double'))`.
- Lines 46-47: Call GRAPE and request analytical gradient; implemented by `[~,~,grad_anl]=grape_phase(guess,spin_system)`.
- Lines 50-51: Left waveform edge; implemented by `wave_forw=guess; wave_forw(1)=wave_forw(1)+h`.
- Lines 62-63: Right waveform edge; implemented by `wave_forw=guess; wave_forw(end)=wave_forw(end)+h`.
- Lines 74-75: Waveform midpoint; implemented by `wave_forw=guess; wave_forw(3)=wave_forw(3)+h`.

### Control flow inferred from the code

- Line 14: `for` loop over `n=1:numel(formalisms)`.
- Line 56: conditional branch on `abs(grad_anl(1)-grad_num)/abs(grad_num)<1e-6`.
- Line 68: conditional branch on `abs(grad_anl(end)-grad_num)/abs(grad_num)<1e-6`.
- Line 80: conditional branch on `abs(grad_anl(3)-grad_num)/abs(grad_num)<1e-6`.

### Key state/data transformations

- Lines 11: computes `formalisms` using `formalisms={'sphten-liouv','zeeman-liouv','zeeman-hilb'}`.
- Lines 17: computes `[spin_system,Sx,Sy,Sz,Lx,Ly,H]` using `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalisms{n})`.
- Lines 20: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 21: computes `control.channels` using `control.channels=[1;1]`.
- Lines 22: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 23: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 24: computes `control.rho_init` using `control.rho_init={ Sx Sy Sz}`.
- Lines 25: computes `control.rho_targ` using `control.rho_targ={-Sz Sy Sx}`.
- Lines 26: computes `control.pwr_levels` using `control.pwr_levels=2*pi*linspace(50e3,70e3,10)`.
- Lines 27: computes `control.method` using `control.method='lbfgs'`.
- Lines 28: computes `control.max_iter` using `control.max_iter=1000`.
- Lines 29: computes `control.plotting` using `control.plotting={}`.
- Lines 30: computes `control.integrator` using `control.integrator='rectangle'`.
- Lines 33: computes `control.distortion` using `control.distortion={@(w)firf(w,[0.9 0.1i]), @(w)spf(w,0.2), @(w)szf(w,0.2), @(w)amp_root(w,2*pi*20e3,4)`.
- Lines 37: computes `control.pulse_dt` using `control.pulse_dt=12.8e-6*ones(1,5)`.
- Lines 38: computes `control.amplitudes` using `control.amplitudes=ones(1,5)`.
- Lines 41: computes `spin_system` using `spin_system=optimcon(spin_system,control)`.
- Lines 44: computes `guess` using `guess=randn(1,5)/3; h=sqrt(eps('double'))`.

## Implementation structure

- Directional derivative test for the phase-modulated GRAPE
- module, with an ensemble of chained linear and non-linear
- instrumental distortions of the waveform.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Define an ensemble of distortion chains
- Set the interval grid
- Spinach housekeeping
- Random phases and finite diff increment
- Call GRAPE and request analytical gradient

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `firf()`, `spf()`, `szf()`, `amp_root()`, `optimcon()`, `eps()`, `grape_phase()`, `squeeze()`, `grad_anl()`, `wave_forw()`, `wave_back()`, `fid_forw()`, `fid_back()`.
