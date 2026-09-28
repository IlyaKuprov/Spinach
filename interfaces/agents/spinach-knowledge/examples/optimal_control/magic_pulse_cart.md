# examples/optimal_control/magic_pulse_cart.m

- Signature: `magic_pulse_cart()`

## Purpose

A template for optimising a broadband ¹³C 90-degree “magic pulse” that tolerates resonance offsets and RF power calibration errors. See [the cited magic-pulse paper](http://dx.doi.org/10.1016/j.jmr.2005.12.010).

At 28.18 T, the pulse is intended to excite approximately 200 ppm (60 kHz) uniformly. To make the worst-case ¹³C–¹H J-coupling (about 200 Hz) negligible, its duration is capped at `1/(100*J) = 50 µs`. The required transfers are `{Lz → Lx, Ly → Ly, Lx → −Lz}`; the anticipated nutation-frequency range across the RF coil is 50–70 kHz. Calculation time: minutes.

## Physical / mathematical content

The example models 100 non-interacting ¹³C spins at equally spaced chemical shifts from −100 to +100 ppm. It uses a spherical-tensor Liouville-space basis with `IK-2` approximation, proximity level 1, and scalar-coupling connectivity. The `Lx`, `Ly`, and `Lz` starting states are normalised before optimisation; their targets are `−Lz`, `Ly`, and `Lx`, respectively.

## Numerical / algorithmic content

- Cartesian RF controls use the `Lx` and `Ly` operators, mapped to the ¹³C channel. The pulse has 40 intervals of 1 µs each, with ten power levels spanning `2π × 50–70 kHz`.
- GRAPE optimisation calls `fmaxnewton` with `@grape_xy` and the `lbfgs` method. The initial guess is a `2 × 40` array of `1/4`; penalties `NS` and `SNS` have weights `0.01` and `10.0`, and the iteration limit is 200. Requested plots are `xy_controls`, `robustness`, and `spectrogram`.
- The optimised profile is scaled by the mean power level and simulated as an XY-shaped pulse using `expv-pwc`. A ¹³C free induction decay is acquired with a 70,000 Hz sweep, 2,048 points, and 16,384-point zero filling; Gaussian apodisation with parameter 10 precedes the Fourier transform. The real spectrum is plotted on an inverted ppm axis.
- For comparison, the script simulates and plots a conventional hard pulse at zero offset, phase `π/2`, power `2π × 60 kHz`, duration `4.2 µs`, rank 3, and `expv` propagation, using the same acquisition and spectral processing.

## Implementation structure

The function sets up the spin system and basis, constructs states and control operators, configures the optimisation, extracts the Cartesian waveform, and compares simulated spectra from the optimised and conventional pulses.

Source contacts: ilya.kuprov@weizmann.ac.il; david.goodwin@inano.au.dk.