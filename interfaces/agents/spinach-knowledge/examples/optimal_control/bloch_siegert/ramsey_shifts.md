# examples/optimal_control/bloch_siegert/ramsey_shifts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/ramsey_shifts.m`
- Signature: `ramsey_shifts()`
- Total lines: 107

## Purpose

Ramsey shifts of other spins under an off-resonant drive. A proton channel drives a three-spin system in which 13C and 15N are far off resonance: each of their Zeeman frequencies shifts by delta_n=w_n*w1n^2/(w_n^2-w_c^2) where w_n is the signed Zeeman frequency of the nucleus, w_c is the signed carrier frequency, and w1n is the drive amplitude scaled by the ratio of the magnetogyric ratios. The negative magnetogyric 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `shift_phases()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Ramsey shifts of other spins under an off-resonant drive. A proton
- channel drives a three-spin system in which 13C and 15N are far off
- resonance: each of their Zeeman frequencies shifts by
- delta_n=w_n*w1n^2/(w_n^2-w_c^2)
- where w_n is the signed Zeeman frequency of the nucleus, w_c is the
- signed carrier frequency, and w1n is the drive amplitude scaled by
- the ratio of the magnetogyric ratios. The negative magnetogyric
- ratio of 15N makes its shift precess the opposite way to 13C. The
- accumulated phases of the off-resonant nuclei are compared with the
- analytic formula, and the quadratic scaling of the shift in the
- drive amplitude and its inverse scaling in the magnet field are
- demonstrated.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `shift_phases()`, `num2str()`, `sign()`, `create()`, `basis()`, `state()`, `operator()`, `true()`, `optimcon()`, `ensemble()`, `bfrq()`, `gam()`.
