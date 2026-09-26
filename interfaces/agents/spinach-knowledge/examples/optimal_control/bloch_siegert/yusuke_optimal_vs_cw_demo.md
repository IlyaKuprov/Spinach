# examples/optimal_control/bloch_siegert/yusuke_optimal_vs_cw_demo.m

- Signature: `yusuke_optimal_vs_cw_demo()`

## Purpose

Bloch-Siegert-aware phase optimisation compared to a simple constant- phase low-power cycle. This is the control-side companion to yusuke_14n_broadening_demo.m: the task is a reduced-model surrogate for low-power offset-tolerant decoupling, formulated as an identity cycle that should preserve magnetisation across offset and B1 distributions. The "non-optimal pulse" is a constant-phase X pulse with the same RF amplitu

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Bloch-Siegert-aware phase optimisation compared to a simple constant-
- phase low-power cycle. This is the control-side companion to
- yusuke_14n_broadening_demo.m: the task is a reduced-model surrogate for
- low-power offset-tolerant decoupling, formulated as an identity cycle
- that should preserve magnetisation across offset and B1 distributions.
- The "non-optimal pulse" is a constant-phase X pulse with the same RF
- amplitude and total duration. It behaves like a simple CW-style cycle:
- acceptable near the design point, but poor once offset and B1 errors
- are included. The "optimal pulse" is a phase-modulated waveform
- optimised with Bloch-Siegert corrections enabled.
- Magnetic field corresponding to 800 MHz 1H
- Single-spin surrogate model
