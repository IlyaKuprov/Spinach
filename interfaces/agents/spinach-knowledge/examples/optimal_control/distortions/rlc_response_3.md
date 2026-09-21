# examples/optimal_control/distortions/rlc_response_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/rlc_response_3.m`
- Signature: `rlc_response_3()`
- Total lines: 102

## Purpose

Probe circuit response effect on the accuracy of the deu- terium pre-phasing pulse designed to set deuterium magne- tisation in a -CD3 group of alanine up for rephasing 100 microseconds after the pulse is finished. The system is assumed to be a powder (100 orientations) with a B1 distribution (from 40 to 60 kHz per channel). Piecewise-linear GRAPE pulse is used. Calculation time: minutes

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Probe circuit response effect on the accuracy of the deu-
- terium pre-phasing pulse designed to set deuterium magne-
- tisation in a -CD3 group of alanine up for rephasing 100
- microseconds after the pulse is finished.
- The system is assumed to be a powder (100 orientations)
- with a B1 distribution (from 40 to 60 kHz per channel).
- Piecewise-linear GRAPE pulse is used.
- Calculation time: minutes
- 600 MHz magnet
- Isotopes
- Alanine CD3 NQI parameters
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `anas2mat()`, `create()`, `basis()`, `drifts()`, `state()`, `operator()`, `false()`, `optimcon()`, `guess()`, `fmaxnewton()`, `pulse_profile()`, `kfigure()`, `restrans()`, `spin()`.
