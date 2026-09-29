# examples/fundamentals/derivative_tests/dirdiff_4_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_4_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_4_rect.m)

- Signature: `dirdiff_4_rect()`

## Purpose

This example checks selected phase-waveform derivatives from the phase-modulated GRAPE routine `grape_phase` against centred finite differences of its reported fidelity, using the rectangle integrator.

## Spin system and controls

For each of `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, `dirdiff_test_system` supplies a test system and the operators `Sx`, `Sy`, `Sz`, `Lx`, `Ly`, and `H`. The control structure sets isotope `13C`, channel map `[1;1]`, drift `H`, controls `Lx` and `Ly`, initial states `{Sx Sy Sz}`, target states `{-Sz Sy Sx}`, and power levels `2*pi*linspace(50e3,70e3,10)`. It selects `method='lbfgs'`, `max_iter=1000`, empty plotting options, `integrator='rectangle'`, `pulse_dt=12.8e-6*ones(1,5)`, and fixed amplitudes `ones(1,5)`, then passes the system and controls through `optimcon`.

The L-BFGS method and iteration limit are control settings here; the example does not run an optimisation loop. It evaluates `grape_phase` directly on a random five-element phase waveform `randn(1,5)/3`. The analytical gradient is taken from the first fidelity component (`grad_anl(:,:,1)`).

## Finite-difference checks

The five-element phase vector is perturbed with `h=sqrt(eps('double'))`. At indices 1, 3, and 5, the script compares the analytic phase gradient with a central difference of the first fidelity component, `(f_1(phi+h*e_k)-f_1(phi-h*e_k))/(2*h)`. Each relative discrepancy must be below `1e-6` or the script raises an error; these are the first, middle, and last phase coordinates.

## Scope

The example checks three coordinates for each constructed formalism; it does not test every phase coordinate or demonstrate optimiser convergence. The source does not declare an explicit unit or phase convention for the waveform values. Its relative-error expression also has no small-denominator guard when a finite-difference estimate is zero or near zero.