# examples/optimal_control/state_transfer_coop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_coop.m`
- Signature: `state_transfer_coop()`
- Total lines: 89

## Purpose

Optimal control pulse optimisation for state-to-state transfer in a quadrupolar 14N spin at a fixed orientation and power level Two pulses are optimised cooperatively, so that the sum of the outcomes only contains the target state, and no impurities. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer
- in a quadrupolar 14N spin at a fixed orientation and power level
- Two pulses are optimised cooperatively, so that the sum of the
- outcomes only contains the target state, and no impurities.
- Calculation time: minutes.
- Magnet field
- Isotopes
- Glycine NQI, random orientation
- Glycine 14N chemical shift
- Basis set
- Run Spinach housekeeping
- Set up and normalise the initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `orientation()`, `carrier()`, `rotframe()`, `operator()`, `optimcon()`, `fmaxnewton()`, `outcome_b()`.
