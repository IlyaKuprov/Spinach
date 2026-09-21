# examples/optimal_control/bloch_siegert/coote_goodcop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/coote_goodcop.m`
- Signature: `coote_goodcop()`
- Total lines: 147

## Purpose

Reproduction of the GOODCOP pulse design logic from Coote et al. with Bloch-Siegert corrections enabled in the optimiser and simulator The pulse enforces contracted-time C-alpha evolution while inverting CO

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `ppm2hz()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Reproduction of the GOODCOP pulse design logic from Coote et al.
- with Bloch-Siegert corrections enabled in the optimiser and simulator
- The pulse enforces contracted-time C-alpha evolution while inverting CO
- Magnetic field corresponding to 800 MHz 1H
- Single-spin carbon model
- Basis set
- Spinach housekeeping
- Relevant operators and states
- Drift Hamiltonian
- Paper parameters for GOODCOP
- Offset grids from the paper
- Convert to offset frequencies

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `operator()`, `state()`, `hamiltonian()`, `assume()`, `ppm2hz()`, `step()`, `ca_hz()`, `true()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `false()`, `eval_hz()`, `bloch_siegert()`.
