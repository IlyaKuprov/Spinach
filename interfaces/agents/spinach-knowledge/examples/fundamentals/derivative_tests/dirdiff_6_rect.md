# examples/fundamentals/derivative_tests/dirdiff_6_rect.m

- Signature: `dirdiff_6_rect()`

## Purpose

Check the analytical Cartesian GRAPE Hessian against finite-differenced gradients.

## Physical / mathematical content

The test checks selected Hessian columns for a five-interval Cartesian waveform in the `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb` formalisms.

## Numerical / algorithmic content

Each numerical column is formed by a centered difference of GRAPE gradients with increment `sqrt(eps('double'))`. The leftmost, middle, and rightmost columns are checked using a relative 1-norm tolerance of `1e-6`, scaled by the corresponding numerical column norm.

## Implementation structure

- Configure the rectangular-integrator, Newton-method GRAPE control system and obtain the analytical Hessian from `grape_xy`.
- Perturb waveform coordinates 1, 3, and 10 in turn, evaluate gradients at both perturbed waveforms, and form the corresponding numerical Hessian columns.
- Compare those columns with the analytical Hessian and raise an error if any comparison exceeds tolerance.
