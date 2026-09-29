# examples/fundamentals/derivative_tests/dirdiff_3_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_3_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_3_rect.m)

- Signature: `dirdiff_3_rect()`

## Purpose

This example checks selected analytical waveform derivatives from the Cartesian GRAPE routine `grape_xy` against centred finite differences of its reported fidelity, using the rectangle integrator.

## Spin system and controls

For each of `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, `dirdiff_test_system` supplies a test system and the operators `Sx`, `Sy`, `Sz`, `Lx`, `Ly`, and `H`. The control structure sets isotope `13C`, channel map `[1;1]`, drift `H`, controls `Lx` and `Ly`, initial states `{Sx Sy Sz}`, target states `{-Sz Sy Sx}`, and power levels `2*pi*linspace(50e3,70e3,10)`. It selects `method='lbfgs'`, `max_iter=1000`, empty plotting options, `integrator='rectangle'`, and `pulse_dt=12.8e-6*ones(1,5)`, then passes the system and controls through `optimcon`.

The L-BFGS method and iteration limit are control settings here; the example does not run an optimisation loop. It evaluates `grape_xy` directly on a random (2×5) waveform `randn(2,5)/3`. The analytical gradient is taken from the first fidelity component (`grad_anl(:,:,1)`).

## Finite-difference checks

With `h=sqrt(eps('double'))`, the script central-differences the first fidelity component at linear control indices 1, 3, and 10. For each index `k`, it compares the analytic gradient with `(f_1(x+h*e_k)-f_1(x-h*e_k))/(2*h)`; the relative discrepancy must be below `1e-6` or it raises an error.

The source calls these indices left edge, midpoint, and right edge. The control array is 2-by-5, so MATLAB linear index 3 is `guess(1,2)`, not the central time column; index 10 is `guess(2,5)`.

## Scope

These are three sampled coordinate-derivative checks for each constructed formalism, not a comparison over all ten waveform coordinates or an optimisation-convergence test. The relative-error expression has no small-denominator guard when the finite-difference estimate is zero or near zero.