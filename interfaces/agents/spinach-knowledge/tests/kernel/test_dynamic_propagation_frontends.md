# tests/kernel/test_dynamic_propagation_frontends.m

- Signature: `result=test_dynamic_propagation_frontends()`

Tests one-spin propagation against finite-dimensional references. `propagator` is compared with `expm`; `step` covers numeric, zero-time, quadrature, callback, and Hilbert-space branches. `evolution` checks final, trajectory, observable, multichannel, refocus, and total outputs. `krylov` checks final, trajectory, observable, multichannel, and refocus outputs. `reduce` checks projector orthogonality, state preservation, and disabled reduction.
