# examples/optimal_control/magic_pulse_phase.m

- Signature: `magic_pulse_phase()`

## Purpose

A template file for the "magic pulse" optimisations. The term refers to a family of broadband NMR pulses that are tolerant to resonance offsets and power calibration errors: Consider a 13C 90-degree excitation pulse in a 28.18 Tesla magnet. The pulse must uniformly excite a bandwidth of around 200 ppm (60 kHz) and must be short enough for the worst-case 13C-1H J-coupling (ca. 200 Hz) to be negligible. The latter requ

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A template file for the "magic pulse" optimisations. The term refers to
- a family of broadband NMR pulses that are tolerant to resonance offsets
- and power calibration errors:
- Consider a 13C 90-degree excitation pulse in a 28.18 Tesla magnet. The
- pulse must uniformly excite a bandwidth of around 200 ppm (60 kHz) and
- must be short enough for the worst-case 13C-1H J-coupling (ca. 200 Hz)
- to be negligible. The latter requirement caps the duration at 1/100*J
- = 50 us. The pulse must accomplish the following transfers: {Lz -> Lx,
- Ly -> Ly, Lx -> -Lz}. A realistically achievable nutation frequency is
- between 50 kHz and 70 kHz across the RF coil.
- Calculation time: minutes.
- Set the magnetic field
