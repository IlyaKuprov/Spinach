# examples/optimal_control/state_transfer_pro.m

- Signature: `state_transfer_pro()`

## Purpose

Optimal control optimisation of a pulse performing magnetisa- tion transfer from H(N) to C(O) in a typical protein backbone spin system (literature data for shifts and couplings) with a range of pulse powers emulating B1 inhomogeneity and a range of offsets to account for imperfect transmitter placement. The waveform is optimized with LBFGS-GRAPE algorithm with po- int-by-point variation and a penalty on the pulse am

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control optimisation of a pulse performing magnetisa-
- tion transfer from H(N) to C(O) in a typical protein backbone
- spin system (literature data for shifts and couplings) with a
- range of pulse powers emulating B1 inhomogeneity and a range
- of offsets to account for imperfect transmitter placement.
- The waveform is optimized with LBFGS-GRAPE algorithm with po-
- int-by-point variation and a penalty on the pulse amplitude.
- Calculation time: hours.
- Magnetic field
- Spin system
- Textbook chemical shifts, ppm
- Scalar couplings, Hz (literature values)
