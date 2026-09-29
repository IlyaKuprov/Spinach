# tests/kernel/test_dynamic_fp_contexts.m

**Source**: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_fp_contexts.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_fp_contexts.m)

## Purpose

Regression test for the compact `imaging()` and `meshflow()` context hand-off paths. It verifies that both contexts assemble finite, correctly sized generators and phantom-derived initial and detection states, and that spatial flow is conserved.

## What the suite checks

- The `imaging` context on a one-dimensional periodic grid and the `meshflow` context on a closed finite-volume mesh must pass the spin basis dimension and the product of spin and spatial dimensions into the pulse sequence. Their assembled generators must have that product dimension, be square and contain finite entries; the constructed initial state and detection coil must cover the active spatial cells.
- A one-dimensional imaging setup must leave the transverse gradient operators empty. Both flow operators must conserve total spatial mass: their column sums vanish, within the test tolerances, for periodic imaging and closed-boundary mesh diffusion.
- The checks exercise the production context hand-off and report structural and conservation assertions. They do not establish a particular image, spectrum or timed numerical benchmark.

## Inputs and outputs

```matlab
result = test_dynamic_fp_contexts()
```

- **Output**: `result` — regression test result structure with explanatory messages, accumulated through `test_close` and `test_true` assertions.
- Takes no inputs.

## References

- `imaging` — Cartesian-grid imaging context exercised by this test.
- `meshflow` — unstructured-mesh flow context exercised by this test.
- `new_test_result`, `test_close`, `test_true`, `test_spin_system`, `state` — test harness and spin-system utilities used by this test.
