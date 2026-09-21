# examples/optimal_control/state_transfer_wf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_wf.m`
- Signature: `state_transfer_wf()`
- Total lines: 92

## Purpose

A transfer of population from the lowermost energy level in a four-spin system to the uppermost energy level using the wave- function space version of the GRAPE algorithm. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A transfer of population from the lowermost energy level in a
- four-spin system to the uppermost energy level using the wave-
- function space version of the GRAPE algorithm.
- Calculation time: minutes.
- Magnetic field
- Isotopes
- Interactions
- Basis set
- Run Spinach housekeeping
- Bottom ground state to start
- Top excited state to finish
- Get the control operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rho_init()`, `rho_targ()`, `operator()`, `hamiltonian()`, `assume()`, `frqoffset()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
