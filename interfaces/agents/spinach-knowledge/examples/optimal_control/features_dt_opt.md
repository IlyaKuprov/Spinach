# examples/optimal_control/features_dt_opt.m

- Signature: `features_dt_opt()`

## Purpose

Optimisation of slice durations in a composite inversion pulse with specified amplitudes, phases, and a constrain- ed overall duration. The initial guess is 270(-x)360(x)90(y)270(-y)360(y)90(x) [Fig. 3] from https://doi.org/10.1016/0022-2364(83)90133-6 --the optimisation demonstrates that a slightly better pulse of the same power and duration exists. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Optimisation of slice durations in a composite inversion
- pulse with specified amplitudes, phases, and a constrain-
- ed overall duration. The initial guess is
- 270(-x)360(x)90(y)270(-y)360(y)90(x) [Fig. 3]
- from https://doi.org/10.1016/0022-2364(83)90133-6 --the
- optimisation demonstrates that a slightly better pulse of
- the same power and duration exists.
- Calculation time: minutes.
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals over the area
- that needs to be affected by the pulse (25 kHz either side)
- Select a basis set -IK-2 keeps complete basis on each
