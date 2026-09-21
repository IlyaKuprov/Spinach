# examples/optimal_control/features_freeze.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_freeze.m`
- Signature: `features_freeze()`
- Total lines: 114

## Purpose

Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment spin system. The start- ing state is Z-magnetisation on 1H, the destination state is Z-magneti- sation on 19F. There are six control channels. A freeze condition is specified -there are two periods in the control sequence that the optimisation is not allowed to touch. The waveform is optimised with 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across a
- scalar coupling in a hydrofluorocarbon fragment spin system. The start-
- ing state is Z-magnetisation on 1H, the destination state is Z-magneti-
- sation on 19F. There are six control channels.
- A freeze condition is specified -there are two periods in the control
- sequence that the optimisation is not allowed to touch.
- The waveform is optimised with the Newton-Raphson GRAPE algorithm desc-
- ribed in
- with point-by-point variation and a penalty on the waveform exceeding
- a user-specified power threshold. The initial guess is a random pulse;
- the optimisation typically achieves a fidelity of 0.999999.
- Calculation time: minutes.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `guess()`, `false()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
