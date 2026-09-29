# examples/optimal_control/state_transfer_wf.m

Source: [examples/optimal_control/state_transfer_wf.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/state_transfer_wf.m)

- Signature: `state_transfer_wf()`

## Design objective

The example formulates population transfer from the lowest to the highest energy level of a four-spin system using GRAPE in wave-function space. Unlike the companion singlet-to-carbon example, the state vectors here are explicitly initialised as the last and first entries of a 16-element vector.

## Spin model and states

The model has two `1H` and two `13C` spins with `sys.magnet=14.1`. Scalar Zeeman entries are `[1.5, 2.0, 30.0, 40.0]`; scalar couplings are 1–2: 7.0, 1–3: 150, 2–4: 150, and 3–4: 50. Units for these interaction entries are not annotated in this source. The basis uses `zeeman-wavef` with approximation `none`. The initial vector has its final component set to one, and the target has its first component set to one. The source calls these variables `rho_init` and `rho_targ`, but they are wave-function-space vectors in this setup.

## GRAPE configuration

The x/y controls act on proton and carbon channels; the drift is offset using transmitter settings `[1050, 5285]` (units not specified at that assignment). The configured power levels are `2*pi*[460, 480, 500, 520, 540]`, over 100 slices of 150 microseconds (15 ms). The `NS` and `SNS` penalties have weights 0.1 and 10, the optimiser method is `goodwin`, and the iteration cap is 100. A random 4-by-100 guess is passed to `fmaxnewton` with `@grape_xy`; the returned pulse is scaled by the mean power level.

## Test calculation

The script propagates the pulse with `shaped_pulse_xy` and computes `real(rho_targ'*rho)` for reporting. This source describes an example and a test path; it does not include a measured numeric fidelity or evidence that the optimisation converged. The source estimates calculation time as minutes.
