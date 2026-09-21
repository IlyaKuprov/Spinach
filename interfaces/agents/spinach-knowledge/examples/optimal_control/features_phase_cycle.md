# examples/optimal_control/features_phase_cycle.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_phase_cycle.m`
- Signature: `features_phase_cycle()`
- Total lines: 111

## Purpose

Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment spin system. The start- ing state is Z-magnetisation on 1H, the destination state is quadrature transverse magnetisation on 19F. A phase cycle is specified: a flip in the phase of the fluorine channel must produce the corresponding flip in the phase of the resulting magne- tisation on 19F. Calculati

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across a
- scalar coupling in a hydrofluorocarbon fragment spin system. The start-
- ing state is Z-magnetisation on 1H, the destination state is quadrature
- transverse magnetisation on 19F.
- A phase cycle is specified: a flip in the phase of the fluorine channel
- must produce the corresponding flip in the phase of the resulting magne-
- tisation on 19F.
- Calculation time: minutes.
- Magnetic field
- Spin system
- Chemical shifts, ppm
- Scalar couplings, Hz (literature values)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `pulse()`, `phased_pulse()`, `report()`, `int2str()`, `mat2cell()`, `shaped_pulse_xy()`, `stateinfo()`.
