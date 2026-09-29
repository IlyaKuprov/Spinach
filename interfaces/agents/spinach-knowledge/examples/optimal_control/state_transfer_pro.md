# examples/optimal_control/state_transfer_pro.m

Source: [examples/optimal_control/state_transfer_pro.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/state_transfer_pro.m)

- Signature: `state_transfer_pro()`

## Design objective

Design a shaped pulse that transfers magnetisation from the amide `1H` (spin 2) to the carbonyl `13C` (spin 5) in a six-spin protein-backbone model. The source describes a robustness ensemble for transmitter placement and RF-amplitude variation, optimised with point-by-point LBFGS-GRAPE and an amplitude penalty.

## Spin model and states

The model uses isotopes `15N, 1H, 13C, 13C, 13C, 15N`, a field setting of `sys.magnet=9.4`, and the source's textbook shift values `[119.79, 8.03, 57.32, 27.71, 177.25, 115.55]` ppm. Its nonzero scalar couplings are 1–3: −11 Hz, 2–3: 140 Hz, 3–4: 35 Hz, 3–5: 55 Hz, 3–6: 7 Hz, and 5–6: −15 Hz. The basis is spherical-tensor Liouville space (`sphten-liouv`) with `IK-0` approximation. The normalised initial and target operators are `Lz` on spins 2 and 5, respectively.

## Robust pulse design

Three channels control `1H`, `13C`, and `15N`, with x/y controls for each. The ensemble contains three offsets per channel on the grid −100, 0, and 100 Hz. The drift Hamiltonian is shifted using source settings `[3214, 10000, -4800]` for the three transmitters; the source does not annotate units on this separate transmitter-placement vector. Five RF power levels span `2*pi*0.9e3` to `2*pi*1.1e3 rad/s`. The waveform has 500 slices of 40 microseconds each (20 ms total), uses the `NS` penalty with weight 0.01, and allows up to 500 LBFGS iterations.

The random 6-by-500 initial guess is seeded with a short y-control segment on the proton channel and a terminal y-control segment on the carbon channel; these are starting-guess features, not a claim about the optimised waveform. `fmaxnewton` is called with `@grape_xy`; the resulting pulse is scaled by the mean power level and propagated with `shaped_pulse_xy`.

## What the example reports

The script computes and reports `real(rho_targ'*rho)` after the test propagation. The source gives no numerical value for the computed overlap or optimisation fidelity. The source estimates calculation time as hours. It cites [DOI 10.1016/j.jmr.2011.07.023](https://doi.org/10.1016/j.jmr.2011.07.023).
