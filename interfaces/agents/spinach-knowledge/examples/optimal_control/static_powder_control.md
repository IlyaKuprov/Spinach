# examples/optimal_control/static_powder_control.m

- Signature: `static_powder_control()`

## Purpose

Optimal control optimisation for a pulse that is designed to set deuterium magnetisation in a -CD3 group of alanine up for perfect rephasing 100 microseconds after the pulse is finished. The system is a powder (100 orientations) with a B1 dist- ribution (from 46 to 54 kHz per channel) and transmitter offset error within 1 kHz of the chemical shift. Goodwin's very efficient version of the GRAPE Hessian al- gorithm is 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Optimal control optimisation for a pulse that is designed
- to set deuterium magnetisation in a -CD3 group of alanine
- up for perfect rephasing 100 microseconds after the pulse
- is finished.
- The system is a powder (100 orientations) with a B1 dist-
- ribution (from 46 to 54 kHz per channel) and transmitter
- offset error within 1 kHz of the chemical shift.
- Goodwin's very efficient version of the GRAPE Hessian al-
- gorithm is used because propagator dimensions are small;
- it yields a sophisticated kind of spin echo.
- Calculation time: minutes
- 600 MHz magnet
