# examples/optimal_control/features_newton.m

- Signature: `features_newton()`

## Purpose

Optimal control pulse optimisation for state-to-state transfer across two scalar couplings in a hydrofluorocarbon fragment spin system. The starting state is Lz on 1H, the destination state is Lz on 19F. There are six control channels, the pulse is designed to be stable with res- pect to proton transmitter offset and pulse nutation frequency drop. The optimisation uses Newton-Raphson GRAPE algorithm described in: wit

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across
- two scalar couplings in a hydrofluorocarbon fragment spin system. The
- starting state is Lz on 1H, the destination state is Lz on 19F. There
- are six control channels, the pulse is designed to be stable with res-
- pect to proton transmitter offset and pulse nutation frequency drop.
- The optimisation uses Newton-Raphson GRAPE algorithm described in:
- with point-by-point variation and a penalty on the waveform exceeding
- a user-specified power threshold. The initial guess is a random pulse.
- Calculation time: minutes.
- Magnetic field
- Spin system
- Chemical shifts, ppm
