# tests/kernel/test_dynamic_trajectory_frontends.m

- Signature: `result=test_dynamic_trajectory_frontends()`

## Purpose

Regression tests for dynamic-trajectory comparison frontends.

## Tests

- Constructs test trajectory data directly from `unit_state()` and `state()` calls; the suite does not propagate dynamics.
- Checks trajectory views for coherence order, per-spin components, and level populations, including explicit time-axis handling.
- Checks `trajsimil()` on the constructed trajectories.

## Outputs

- result -regression test result with explanatory messages
- The test exercises trajan() plotting branches and trajsimil() scoring
- branches on a compact two-spin spherical-tensor trajectory.
