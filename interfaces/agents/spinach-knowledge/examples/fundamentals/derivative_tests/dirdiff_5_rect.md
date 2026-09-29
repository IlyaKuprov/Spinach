# examples/fundamentals/derivative_tests/dirdiff_5_rect.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_5_rect.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_5_rect.m)

- Signature: `dirdiff_5_rect()`

## Question tested

This is an internal-consistency comparison between the Newton and Goodwin GRAPE Hessians, not a finite-difference test. It checks both phase-modulated and Cartesian controls with the rectangular integrator.

## Setup

For each of `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, `dirdiff_test_system` supplies the spin system, states, operators, and drift. The controls use isotope `13C`, channel map `[1;1]`, drift `H`, controls `Lx` and `Ly`, initial states {Sx,Sy,Sz}, targets {-Sz,Sy,Sx}, and power levels `2*pi*linspace(50e3,70e3,10)`. The configuration also sets `max_iter=1000` and an empty plotting list. The rectangular grid is `12.8e-6*ones(1,5)`, and the amplitude vector has five entries. A phase guess `randn(1,5)/3` is used for `grape_phase`; a separate Cartesian guess `randn(2,5)/3` is used for `grape_xy`.

For each modulation, the script obtains one Hessian with `control.method='newton'` and another with `control.method='goodwin'`, using the same waveform guess for both methods within that modulation.

## Comparison and observable

For each pair, the Hessian arrays are flattened and their one-norm difference is compared with the one-norm of the Newton Hessian. The script raises an error when `norm(H_newton(:)-H_goodwin(:),1) > 1e-6*norm(H_newton(:),1)`. Otherwise it prints a formalism-specific consistency message. This comparison is performed once for the phase Hessians and once for the Cartesian Hessians.

## Scope

The check establishes only whether these two Hessian outputs agree within the coded relative one-norm threshold for these two random guesses and the configured control system. It does not compare either Hessian with an independent finite-difference reference.