# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1.m

- Signature: `xix_q_con_time_ensemble_b1()`

## Purpose

Calculates steady-state proton polarisation versus XiX contact time while averaging over an ensemble of electron microwave-field (Rabi-frequency) values.

## Model and ensemble

The source models a trityl electron–proton pair at 1.2142 T and 80 K, with 3.5 nm separation and orientation-dependent proton relaxation from `r1n_dnp`. It uses T1/T2 relaxation, diagonal relaxation terms and the `dibari` equilibrium, with the full spherical-tensor Liouville basis and no basis approximation. Five Gauss–Legendre nodes span microwave-field values of 10–20 MHz; their weights are used to average the calculated proton signal.

## XiX scan and output

The experiment detects proton `Lz` on an 800-point two-angle spherical powder grid. It scans 1–64 XiX loops using 48 ns pulses (the second pulse has inverted phase); the contact time is twice the loop count times the pulse duration. For each field node and loop count, the steady state is computed with `powder(...,@xixdnp_steady,...,'esr')`. The source sets 153 μs shot spacing minus the total pulse duration, and uses −13 MHz added shift and +61 MHz electron offset. The field-weighted real proton expectation value is plotted against contact time and saved as `xix_q_con_time_ensemble_b1.fig`.
