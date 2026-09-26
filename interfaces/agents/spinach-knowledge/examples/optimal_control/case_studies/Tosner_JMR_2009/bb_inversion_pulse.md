# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m

- Signature: `bb_inversion_pulse()`

## Purpose

Broadband inversion pulse design for liquid-state NMR. Reprodu- ces, using Spinach, the second example from: A single proton is considered in the rotating frame with multiple transmitter offsets (or chemical shifts). The goal is to design a 600 µs broadband inversion pulse (1 µs slices) that performs: I_z → -I_z uniformly over a frequency offset range of ±50 kHz; controls are Cartesian (Lx, Ly) operators in the rotat

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Broadband inversion pulse design for liquid-state NMR. Reprodu-
- ces, using Spinach, the second example from:
- A single proton is considered in the rotating frame with multiple
- transmitter offsets (or chemical shifts). The goal is to design a
- 600 µs broadband inversion pulse (1 µs slices) that performs:
- I_z → -I_z
- uniformly over a frequency offset range of ±50 kHz; controls are
- Cartesian (Lx, Ly) operators in the rotating frame.
- Magnetic field (Tesla)
- Chemical shift (ppm)
- Basis set
- Spinach housekeeping
