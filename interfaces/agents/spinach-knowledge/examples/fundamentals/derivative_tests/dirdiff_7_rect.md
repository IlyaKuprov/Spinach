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
