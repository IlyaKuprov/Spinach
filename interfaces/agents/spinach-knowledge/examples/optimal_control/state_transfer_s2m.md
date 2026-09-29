# examples/optimal_control/state_transfer_s2m.m

Source: [examples/optimal_control/state_transfer_s2m.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/state_transfer_s2m.m)

- Signature: `state_transfer_s2m()`

## Design objective

This example designs a pulse to transfer coherence from a two-proton singlet to a nearby carbon, in a setting described as typical of parahydrogenation. Its source identifies LBFGS-GRAPE and cites [DOI 10.1016/j.jmr.2011.07.023](https://doi.org/10.1016/j.jmr.2011.07.023).

## Spin model and transfer states

The four-spin system contains two `1H` and two `13C` spins, with `sys.magnet=14.1`. The scalar Zeeman entries are `[1.5, 2.0, 30.0, 40.0]`; the source does not label their units. Nonzero scalar couplings are 1–2: 7.0, 1–3: 150, 2–4: 150, and 3–4: 50 (the source gives no unit annotation beside these entries). The basis is `sphten-liouv` with approximation `none`. The normalised starting state is the singlet on spins 1 and 2; the normalised target is `Lz` on spin 4.

## Pulse and ensemble settings

The controls are x/y operators on the proton and carbon channels. The drift Hamiltonian is shifted with transmitter settings `[1050, 5285]` (units are not annotated at that assignment). Five power levels are specified as `2*pi*[460, 480, 500, 520, 540]`; the waveform has 100 slices of 150 microseconds each, a 15 ms design duration. LBFGS is configured for at most 100 iterations, with `NS` and `SNS` penalties weighted 0.1 and 10. A random 4-by-100 guess is used, then `fmaxnewton` optimises with `@grape_xy`; the pulse is scaled by the mean configured power before propagation.

## Reported quantity and evidence boundary

The test simulation computes `real(rho_targ'*rho)` and sends the value to Spinach's report function. The source preamble reports a 50% terminal fidelity as a benchmark; the file embeds no numerical output from its test propagation. The source estimates minutes for calculation time.
