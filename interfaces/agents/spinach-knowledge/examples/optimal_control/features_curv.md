# examples/optimal_control/features_curv.m

- Signature: `features_curv()`

## Purpose

A transfer of coherence from longitudinal magnetization into a two-spin singlet state with a distribution of B1 powers. An ensemble of ten spin systems with different power levels is simultaneously driven to optimal fidelity, which in this case is 1/sqrt(2) = 0.7071 Curvilinear GRAPE interface is used -the user specifies the definition of the curvilinear coordinates and the Jacobian. In this case, the coor- dinates a

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A transfer of coherence from longitudinal magnetization into a two-spin
- singlet state with a distribution of B1 powers. An ensemble of ten spin
- systems with different power levels is simultaneously driven to optimal
- fidelity, which in this case is 1/sqrt(2) = 0.7071
- Curvilinear GRAPE interface is used -the user specifies the definition
- of the curvilinear coordinates and the Jacobian. In this case, the coor-
- dinates are phase-amplitude.
- Calculation time: minutes.
- Magnetic field
- Isotopes
- Interactions
- Basis set
