# examples/optimal_control/pattern_pulse_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/pattern_pulse_1.m`
- Signature: `pattern_pulse_1()`
- Total lines: 103

## Purpose

Nutation frequency selective excitation described in Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified nutation frequency intervals have magnetisation arriving into us- er specified states. The pulse is phase-modulated. Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Nutation frequency selective excitation described in Glaser group
- paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified
- nutation frequency intervals have magnetisation arriving into us-
- er specified states. The pulse is phase-modulated.
- Calculation time: minutes.
- Magnetic field
- Single carbon spin
- Transmitter is on resonance
- No approximations
- Run Spinach housekeeping
- Get pertinent spin states
- Get pertinent control operators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `figure()`, `kxlabel()`, `klegend()`, `ylim()`, `num2cell()`, `optimcon()`, `fmaxnewton()`, `polar2cartesian()`, `shaped_pulse_xy()`.
