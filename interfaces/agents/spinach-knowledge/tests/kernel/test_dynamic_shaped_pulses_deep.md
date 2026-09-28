# tests/kernel/test_dynamic_shaped_pulses_deep.m

- Signature: `result=test_dynamic_shaped_pulses_deep()`

## Purpose

Tests dynamic shaped-pulse propagation paths. Syntax: result=test_dynamic_shaped_pulses_deep()

## Physical / mathematical content

- Tests shaped-pulse propagators against exact references for a constant RF generator in a one-spin Liouville-space system. The checks compare final states and trajectories.
## Numerical / algorithmic content

- Exercises `shaped_pulse_xy` with constant piecewise-constant and piecewise-linear controls using `expv`, `expm`, and evolution-product paths. It also exercises `shaped_pulse_af` with constant amplitude and zero frequency offset using `expv`, `expm`, and `evolution`.
## Outputs

- `result` — regression test result with explanatory messages.
## Implementation structure

- Build a one-proton Liouville-space spin system and obtain its controls.
- Define a constant Cartesian RF generator over two slices; compare `shaped_pulse_xy` final states and trajectories for `expv-pwc` and `expv-pwl` with the exact constant-generator reference.
- Test the corresponding `expm-pwc`, `expm-pwl`, `evol-pwc`, and `evol-pwl` paths against that reference.
- Define a constant amplitude-frequency pulse over three slices; compare the `shaped_pulse_af` `expv`, `expm`, and `evolution` final states and trajectories with the exact reference, and compare the returned propagator for the `expm` path.