# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m

- Signature: `xix_w_pulse_dur_ensemble_b1()`

## Purpose

2D parameter scan of XiX DNP in the steady state with electron Rabi frequency ensemble. Calculation time: hours.

## Physical / mathematical content

- Models an electron–proton pair at 3.4 T and 80 K, separated by 3.5 Å, with specified electron g-tensor and proton Zeeman shifts. Uses `t1_t2` relaxation, including an orientation-dependent proton R1 rate, and detects the proton Lz expectation value. The XiX sequence uses a phase-inverted second pulse.

## Numerical / algorithmic content

- Samples five electron B1 values from 10 to 20 MHz with Gaussian-Legendre weights. For each B1 value, a `parfor` loop scans 200 pulse durations from 2 to 21 ns; each steady-state calculation uses ESR powder averaging over the specified spherical grid and 101 electron offsets from −230 to 205 MHz. The results are weighted over the B1 ensemble before plotting.

## Implementation structure

- Sets the W-band field, electron and proton interactions, temperature, coordinates, spherical-tensor Liouville-space basis, and propagator tolerance.
- Constructs the spin system and proton detection operator, then sets the XiX experiment parameters, including pulse phase, offset grid, loop count, and shot spacing.
- Calls `powder` with `@xixdnp_steady` for each B1 value and pulse duration; plots the real proton expectation value against microwave offset and pulse duration, then saves `xix_w_pulse_dur_ensemble_b1.fig`.
