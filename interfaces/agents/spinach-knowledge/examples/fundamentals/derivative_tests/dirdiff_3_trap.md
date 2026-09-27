# examples/fundamentals/derivative_tests/dirdiff_3_trap.m

- Signature: `dirdiff_3_trap()`

## Purpose

Check directional derivatives of the Cartesian GRAPE module with the trapezium integrator.

## Physical / mathematical content

The test checks how the GRAPE fidelity changes with selected Cartesian waveform samples. Spin systems are constructed in the `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb` formalisms.

## Numerical / algorithmic content

For the left edge, midpoint, and right edge of a random two-channel waveform, the analytical gradient from `grape_xy` is compared with a centered finite difference using `sqrt(eps('double'))). Each relative discrepancy must be below `1e-6`.

## Implementation structure

- Configure Cartesian controls with the trapezium integrator, L-BFGS method, and `12.8e-6` s pulse intervals.
- Obtain the analytical gradient for a random `2×5` waveform.
- Perturb waveform entries 1, 3, and 10 in both directions and compare the resulting finite-difference gradients with the corresponding analytical entries.
- Raise an error for any failed edge or midpoint check.
