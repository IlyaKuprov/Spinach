# examples/optimal_control/features_bss.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_bss.m`
- Signature: `features_bss()`
- Total lines: 90

## Purpose

Optimal control pulse optimisation with Bloch-Siegert shift corrections switched on. A single proton with a Larmor frequency of 1 MHz is driven at a significant fraction of its Larmor frequency --a regime where the counter-rotating component of the control field shifts the resonance ap- preciably. A 90-degree pulse is optimised using LBFGS-GRAPE algorithm. For comparison, the same pulse is also optimised with the cor

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control pulse optimisation with Bloch-Siegert shift corrections
- switched on. A single proton with a Larmor frequency of 1 MHz is driven
- at a significant fraction of its Larmor frequency --a regime where the
- counter-rotating component of the control field shifts the resonance ap-
- preciably. A 90-degree pulse is optimised using LBFGS-GRAPE algorithm.
- For comparison, the same pulse is also optimised with the corrections
- switched off and then evaluated in the corrected model.
- Calculation time: minutes.
- Larmor frequency of 1 MHz
- Spin system
- Chemical shifts, ppm
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `state()`, `operator()`, `true()`, `optimcon()`, `fmaxnewton()`, `ensemble()`, `false()`, `report()`, `num2str()`.
