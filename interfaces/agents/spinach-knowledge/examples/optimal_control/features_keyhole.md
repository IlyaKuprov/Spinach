# examples/optimal_control/features_keyhole.m

- Signature: `features_keyhole()`
- Source: [examples/optimal_control/features_keyhole.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_keyhole.m)

## Purpose and model

The script designs a simulated transfer from 1H longitudinal magnetisation to 19F in a 9.4 T 1H-13C-19F model. Chemical shifts are zero and the scalar couplings are 140 Hz (1H-13C) and -160 Hz (13C-19F). The model is specified in code; no experimental data are imported.

## Keyhole and controls

Six x/y controls address the three nuclei. The 50 intervals are each 200 microseconds (10 ms total), and five RF-power levels span `2*pi*800` to `2*pi*1200` rad/s. At sample 20, the keyhole callback evaluates `correlation(spin_system,rho,2,'all')`, imposing an intermediate all-spin second-order-correlation condition while the final target remains the 19F longitudinal state. The source comment describes the control sequence as piecewise-linear; the code supplies 50 interval samples and does not specify a separate interpolation object.

The optimisation uses the SNS penalty (weight 100), `lbfgs`, and a 200-iteration limit, with `fmaxnewton` and `@grape_xy`. It starts from a random 6-by-50 waveform and enables correlation-order, per-spin, x/y-control, and spectrogram plots.

## Output and limits

The optimised waveform is scaled by the mean power, propagated with `shaped_pulse_xy` using `expv-pwc`, and assessed by the real final target-state overlap. This check is one propagation, not a reported power-ensemble table. The result is a numerical model simulation; the example does not establish experimental or hardware performance.
