# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m`
- Signature: `bb_refocusing_pulse()`
- Total lines: 97

## Purpose

Spinach implementation of the broadband refocusing example from GRAPE is used to design a 200 µs broadband x-phase π pulse: {Sx -> Sx, Sy -> -Sy, Sz -> -Sz} over an offset range of ±12.5 kHz.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Spinach implementation of the broadband refocusing example from
- GRAPE is used to design a 200 µs broadband x-phase π pulse:
- {Sx -> Sx, Sy -> -Sy, Sz -> -Sz}
- over an offset range of ±12.5 kHz.
- Magnetic field (Tesla)
- Chemical shift (ppm)
- Basis set
- Spinach housekeeping
- Normalised Cartesian basis states
- RF controls and offset operator
- Drift Hamiltonian
- Control data structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `xy_profile()`, `offs_hz()`, `shaped_pulse_xy()`, `fidelities()`, `kfigure()`, `kxlabel()`, `kylabel()`, `xlim()`.
