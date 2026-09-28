# examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1_r.m

- Signature: `xix_w_field_profile_ensemble_b1_r()`

## Purpose

Simulation of a steady-state XiX DNP field profile averaged over electron–proton distance and electron Rabi frequency. Calculation time: minutes.

## Physical / mathematical content

- Models an electron–proton pair at 3.4 T and 80 K, with a trityl electron g-tensor and a proton chemical-shift guess. Distance- and orientation-dependent proton relaxation is supplied by `r1n_dnp`; the detected observable is proton `Lz`.
- Computes the XiX DNP response across microwave resonance offsets from −300 to 300 MHz, then averages over the B1 distribution and the distance distribution. The distance average includes the radial Jacobian factor r².

## Numerical / algorithmic content

- Uses a full `sphten-liouv` basis, a spherical powder grid, and `powder(...,@xixdnp_steady,...,'esr')` to calculate the steady-state response for each distance and electron nutation frequency.
- Samples three distances from 3.5 to 20 Å and five B1 frequencies from 10 to 20 MHz using Gauss–Legendre points and weights. The field profile contains 201 microwave-offset points.

## Implementation structure

- Sets the magnet field, electron and proton interactions, temperature, basis, propagator tolerance, and simulation options.
- For each distance, sets Cartesian spin coordinates and relaxation rates, creates the spin system, and configures proton detection and XiX experiment parameters, including 18 ns pulses, 10 blocks, and a phase-inverted second pulse.
- For each B1 value, runs the powder-averaged steady-state simulation; then integrates over B1 and distance, plots the real proton `Lz` expectation value against microwave offset, and saves `xix_w_field_profile_ensemble_b1_r.fig`.
