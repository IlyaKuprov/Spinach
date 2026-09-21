# examples/optimal_control/bloch_siegert/bloch_siegert_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/bloch_siegert_a.m`
- Signature: `bloch_siegert_a()`
- Total lines: 106

## Purpose

Bloch-Siegert shift compensation functionality demo. The script optimises a 90-degree pulse (Lz -> Lx) for a sing- le spin on resonance. As the control power is increased, Bloch-Siegert shift starts to reduce the fidelity unless it is correctly accounted for. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Bloch-Siegert shift compensation functionality demo. The
- script optimises a 90-degree pulse (Lz -> Lx) for a sing-
- le spin on resonance. As the control power is increased,
- Bloch-Siegert shift starts to reduce the fidelity unless
- it is correctly accounted for.
- Calculation time: minutes.
- Set magnetic field
- Set isotopes
- Set interactions
- Set basis
- Run Spinach housekeeping
- Build and normalise the initial state (Lz)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `spin()`, `pwr_list()`, `false()`, `optimcon()`, `true()`, `fmaxnewton()`, `fid_a()`, `ensemble()`, `fid_b()`, `kfigure()`.
