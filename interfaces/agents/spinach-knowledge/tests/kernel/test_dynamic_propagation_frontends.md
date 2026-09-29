# tests/kernel/test_dynamic_propagation_frontends.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_propagation_frontends.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_propagation_frontends.m)

## Purpose

Regression test for the dynamic propagation front-end kernels of Spinach. The test exercises `propagator()`, `step()`, `evolution()`, `krylov()`, and `reduce()` against direct finite-dimensional propagation references on tiny spin systems.

## What the suite checks

- `propagator` and `step` must agree with direct small-matrix propagation, including the Taylor, numeric Liouville, time-independent and function-handle routes. A zero-duration step preserves the state; the Hilbert density-matrix route matches its commutator reference.
- `evolution` and `krylov` must match explicitly propagated references across their final-state, full-trajectory, single- and multichannel observable, refocusing and total-output modes. These checks guard output meaning and shape as well as numerical agreement on the test systems.
- `reduce` must return orthonormal projectors whose retained subspaces reconstruct the input state; with trajectory-level reduction disabled it returns the identity projector. The suite compares against direct finite-dimensional references, not a full-scale performance run.

## Inputs and outputs

```matlab
result = test_dynamic_propagation_frontends()
```

- **Output:** `result` — regression test result structure with explanatory messages for each checked branch.
- **Input:** none.

## References

- Spinach dynamic propagation kernels: `propagator()`, `step()`, `evolution()`, `krylov()`, `reduce()`.
- [Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach)
