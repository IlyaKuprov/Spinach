# examples/optimal_control/features_dt_var.m

- Signature: `features_dt_var()`

## Purpose

Demonstrates optimal-control pulse design with nonuniform time slices and RF-power robustness. It optimises a shaped pulse for the stated state-transfer task and checks the result by propagating the pulse.

## Physical / mathematical content

The spin system contains 1H, 13C, and 19F channels, with x- and y-phase controls on each nucleus. The objective includes a state-norm (SNS) penalty and is evaluated over five RF-power levels.

## Numerical / algorithmic content

The 50 pulse slices have nonuniform durations. The source uses the lbfgs method with grape_xy gradients through fmaxnewton, and then validates the optimised waveform with shaped_pulse_xy.

## Implementation structure

The code constructs and normalises the initial and target states, generates an initial control guess, and optimises six x/y controls over the five-point RF-power ensemble with an SNS penalty. It then rescales and propagates the waveform using shaped_pulse_xy and reports the fidelity. The pulse has 50 nonuniform time slices.
