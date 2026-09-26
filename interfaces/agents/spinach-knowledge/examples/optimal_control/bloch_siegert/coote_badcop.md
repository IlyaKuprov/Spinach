# examples/optimal_control/bloch_siegert/coote_badcop.m

- Signature: `coote_badcop()`

## Purpose

Reproduction of BADCOP-style selective decoupling logic from Coote et al. with Bloch-Siegert corrections enabled in the optimiser and simulator BADCOP1, BADCOP2, and BADCOP3 are designed and validated

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Reproduction of BADCOP-style selective decoupling logic from Coote et al.
- with Bloch-Siegert corrections enabled in the optimiser and simulator
- BADCOP1, BADCOP2, and BADCOP3 are designed and validated
- Magnetic field corresponding to 800 MHz 1H
- Single-spin carbon model
- Basis set
- Spinach housekeeping
- Relevant operators and states
- Drift Hamiltonian
- Shared paper parameters
- Build all three variants from Table 1
- Design and evaluate each variant
