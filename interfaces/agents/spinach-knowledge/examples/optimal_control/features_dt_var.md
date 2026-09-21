# examples/optimal_control/features_dt_var.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_dt_var.m`
- Signature: `features_dt_var()`
- Total lines: 94

## Purpose

Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment spin system. The start- ing state is Z-magnetisation on 1H, the destination state is Z-magneti- sation on 19F. There are six control channels, the time slice duration is not uniform, and there is a distribution over control powers. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across a
- scalar coupling in a hydrofluorocarbon fragment spin system. The start-
- ing state is Z-magnetisation on 1H, the destination state is Z-magneti-
- sation on 19F. There are six control channels, the time slice duration
- is not uniform, and there is a distribution over control powers.
- Calculation time: minutes.
- Magnetic field
- Spin system
- Chemical shifts, ppm
- Scalar couplings, Hz (literature values)
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
