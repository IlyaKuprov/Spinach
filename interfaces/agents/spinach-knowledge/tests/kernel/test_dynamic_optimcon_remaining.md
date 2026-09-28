# tests/kernel/test_dynamic_optimcon_remaining.m

- Signature: `result=test_dynamic_optimcon_remaining()`

## Purpose

Regression-tests remaining dynamic optimal-control helpers using small deterministic fixtures. Returns a test result with explanatory messages.

## Test coverage

- Checks two-channel waveform distortions: identity mapping, non-orthogonal channel mixing at 60 degrees, causal complex FIR filtering, single-pole and single-zero filtering, and phase-preserving tanh and root-sigmoid amplitude compression. Compares returned Jacobians with centred finite-difference directional derivatives.
- Tests causal FIR kernel estimation using `backslash`, `pinv`, `svd`, and Tikhonov-regularised solver paths, plus `same` alignment.
- Checks BFGS Hessian updates and bad-curvature safeguards, BFGS history reconstruction, LBFGS inverse-Hessian action, Hessian ordering, and regularisation of an indefinite Hessian to positive definiteness.
- Tests frequency-amplitude-phase-time conversion on a supplied grid; instantaneous-frequency recovery for a quadratic-phase signal; and masking of finite-difference stencils containing weak or zero-magnitude samples.
- Checks drift extraction from a two-member context ensemble and trapezium-product auxiliary matrices, including left and right directional-derivative blocks and mixed-derivative matrix sizes. Exercises Cartesian-control plotting offscreen and verifies plotted instantaneous frequency for a linear-phase signal.
- Compares ensemble GRAPE fidelity and gradient with the Cartesian wrapper, checks identity curvilinear coordinates, and compares phase-control gradients with centred finite differences.
- Checks that `fmaxnewton` with zero maximum iterations returns its initial point without iterations or derivative calls. Tests finite Liouville-space GRAPE fidelity, gradient, and Hessian values, gradient agreement with centred finite differences, and Hessian symmetry for a real fidelity. Checks TGRAPE duration gradients against centred finite differences and verifies finite cooperative-GRAPE outputs and gradient shape.

## Implementation structure

The test ensures a parallel pool is available, then runs independent distortion, quasi-Newton, waveform-utility, and GRAPE-family check groups. Local helper functions construct reference waveforms, small control-system fixtures, and finite-difference comparisons.