# examples/optimal_control/state_transfer_m2s.m

- Signature: `state_transfer_m2s()`
- Source: [`examples/optimal_control/state_transfer_m2s.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/state_transfer_m2s.m)

## Purpose

Design a transfer from longitudinal magnetisation to a two-spin singlet-related target in allyl pyruvate under transmitter-offset and pulse-power distributions. The source says the `XX` and `YY` components are omitted from the target because they dephase rapidly in this system, and estimates the calculation time as many hours. These are the example's design assumptions and estimate, not results of a run here.

## Spin model and target

The script loads the allyl-pyruvate proton system, removes the methyl group by retaining the first five spins, and sets `sys.magnet=11.7464` (the source labels the field as 500.13 MHz). It uses the full spherical-tensor Liouville-space basis. The initial state is longitudinal magnetisation on all spins; the target operator is the negative product of the `Lz` operators for spins labelled `Hb` and `Hc`. The source constructs an NMR drift, applies a 2670 transmitter offset, and uses proton `Lx` and `Ly` controls.

## Robust pulse design

The design samples five offsets from −10 to 10 and five configured power levels, `2*pi*[480 490 500 510 520]`; the source does not attach units to those ensemble values. It uses 300 slices of 1 ms each (0.3 s total), penalties `NS` and `SNS` with weights 1 and 10, L-BFGS, and a 3000-iteration termination limit. The initial 2-by-300 Cartesian pulse is built from a 300 Hz cosine and a constant second component, scaled by 0.05. Optimisation is requested through `fmaxnewton(spin_system,@grape_xy,pulse)`; the plotting configuration includes correlation/coherence order, control components, per-spin controls, amplitude, and spectrogram views.

## Scripted pulse and spectrum comparison

The example first simulates and plots a pulse-acquire spectrum from the initial state, using `parameters.sweep=1000`, 2048 points, 4096-point zero filling, and Hz axis units. It then scales the optimised pulse by the mean configured power level, propagates it with `shaped_pulse_xy`, applies homospoil, and computes `real(rho_targ'*rho)` as a fidelity diagnostic. A second pulse-acquire spectrum is generated from that propagated state for comparison. These are operations encoded in the source; no numerical fidelity or successful run is asserted here.
