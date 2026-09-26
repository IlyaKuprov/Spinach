# examples/optimal_control/features_multitarget.m

- Signature: `features_multitarget()`

## Purpose

An example of multi-target optimal control pulse design in the context of singlet state NMR spectroscopy. A pulse is designed that moves TT (carbon-triplet, proton-triplet) into SS (carbon-singlet, proton-singlet) and TS (carbon-triplet, proton-singlet) into ST (carbon-singlet, proton- triplet) simultaneously. The system is assumed to have a distribution in one of the J-couplings. Calculation time: hours.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- An example of multi-target optimal control pulse design in the context
- of singlet state NMR spectroscopy. A pulse is designed that moves TT
- (carbon-triplet, proton-triplet) into SS (carbon-singlet, proton-singlet)
- and TS (carbon-triplet, proton-singlet) into ST (carbon-singlet, proton-
- triplet) simultaneously. The system is assumed to have a distribution in
- one of the J-couplings.
- Calculation time: hours.
- Magnetic field
- Isotopes
- Interactions
- Basis set
- Run Spinach housekeeping
