# examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m`
- Signature: `coherence_transfer()`
- Total lines: 79

## Purpose

The first optimal control example from A heteronuclear two-spin system (1H–13C) with an scalar and both nuclei set on resonance; the goal is to trans- fer transverse magnetisation from proton to carbon: Hx → Cx over a fixed evolution period T = 1/J.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- The first optimal control example from
- A heteronuclear two-spin system (1H–13C) with an scalar
- and both nuclei set on resonance; the goal is to trans-
- fer transverse magnetisation from proton to carbon:
- Hx → Cx
- over a fixed evolution period T = 1/J.
- Magnetic field, Tesla
- Chemical shifts, ppm
- Scalar coupling, Hz
- Basis set
- Spinach housekeeping
- Initial state: Lx on proton (spin 1)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`.
