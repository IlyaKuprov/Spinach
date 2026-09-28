# examples/fundamentals/derivative_tests/dirdiff_4_rect.m

- Signature: `dirdiff_4_rect()`

## Purpose

Check directional derivatives of the phase-modulated GRAPE module with the rectangles integrator.

## Physical / mathematical content

The test checks how the GRAPE fidelity changes with selected phase-waveform samples. Spin systems are constructed in the `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb` formalisms.

## Numerical / algorithmic content

For the left edge, midpoint, and right edge of a random five-sample phase waveform, the analytical gradient from `grape_phase` is compared with a centered finite difference using `sqrt(eps('double'))`. Each relative discrepancy must be below `1e-6`.

## Implementation structure

- Configure phase-modulated controls with the rectangles integrator, L-BFGS method, unit amplitudes, and `12.8e-6` s pulse intervals.
- Obtain the analytical gradient for a random five-element phase waveform.
- Perturb waveform entries 1, 3, and 5 in both directions and compare the finite-difference gradients with the corresponding analytical entries.
- Raise an error for any failed edge or midpoint check.
