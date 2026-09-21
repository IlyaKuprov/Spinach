# examples/optimal_control/features_diss_drift.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_diss_drift.m`
- Signature: `features_diss_drift()`
- Total lines: 135

## Purpose

Optimal control optimisation of a pulse performing magnetisa- tion transfer from H(N) to C(O) in a typical protein backbone spin system (literature data for shifts and couplings) with a range of pulse powers emulating B1 inhomogeneity and a range of offsets to account for imperfect transmitter placement. The dynamics includes dissipative terms in the drift generator: C(O) and N(H) are set to have rapid transverse rel

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
- The dynamics includes dissipative terms in the drift generator:
- C(O) and N(H) are set to have rapid transverse relaxation. Four-
- spin correlation approximation is used, wherein five-spin and
- higher correlations are dropped from the basis set.
- The waveform is optimized with LBFGS-GRAPE algorithm with point-
- by-point variation and a penalty on the waveform amplitude.
- Calculation time: hours.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `assume()`, `hamiltonian()`, `relaxation()`, `frqoffset()`, `optimcon()`, `guess()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
