# examples/fundamentals/derivative_tests/dirdiff_5_rect.m

- Signature: `dirdiff_5_rect()`

## Purpose

Check internal consistency between the Newton and Goodwin GRAPE Hessians.

## Physical / mathematical content

The test compares Hessians for both phase-modulated and Cartesian controls in the `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb` formalisms.

## Numerical / algorithmic content

With the rectangles integrator, the Newton and Goodwin Hessian arrays are compared using their relative 1-norm difference. Each comparison must be no greater than `1e-6` times the Newton Hessian's 1-norm.

## Implementation structure

- Configure the GRAPE control system and generate an initial phase waveform and an initial Cartesian waveform.
- For each waveform, obtain the Hessian from `grape_phase` or `grape_xy` with `control.method` set first to `newton`, then to `goodwin`.
- Compare the two Hessians for each modulation and fail if either consistency check exceeds tolerance.
