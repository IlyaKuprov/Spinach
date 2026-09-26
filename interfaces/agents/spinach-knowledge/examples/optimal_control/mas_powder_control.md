# examples/optimal_control/mas_powder_control.m

- Signature: `mas_powder_control()`

## Purpose

Optimal control pulse starting with Lz and populating the Ly state on 87Rb in a quadrupolar rubidium system under magic angle spinning. A phase-modulated pulse is produced. Calculation time: hours.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Optimal control pulse starting with Lz and populating the
- Ly state on 87Rb in a quadrupolar rubidium system under
- magic angle spinning. A phase-modulated pulse is produced.
- Calculation time: hours.
- System specification
- Quadrupolar coupling
- Basis set and formalism
- Spinach housekeeping
- MAS experiment parameters
- Drift Liouvillians and classical subspace dimension for the ensemble
- Initial state -Lz
- Target state -Ly
