# examples/dnp_sol/steady_state/xix_w_field_profile_single.m

- Signature: `xix_w_field_profile_single()`

## Purpose

Simulate a steady-state XiX DNP field profile for a single electron–proton spin system without ensemble averaging. The source estimates a calculation time of seconds.

## Physical / mathematical content

- The system contains an electron and a proton at 3.4 T and 80 K, with trityl electron g-tensor and estimated proton Zeeman parameters. Their Cartesian coordinates place them 3.5 units apart in the supplied coordinate system. The model uses `t1_t2` relaxation, including a distance- and orientation-dependent proton R1 rate calculated by `r1n_dnp`.

## Numerical / algorithmic content

- The calculation uses an unrestricted `sphten-liouv` basis and diagonal relaxation. It detects proton `Lz` while sweeping 201 electron offsets from −300 to 300 MHz. The experiment specifies a spherical grid, ten XiX blocks, and an inverted second-pulse phase; `powder` runs `xixdnp_steady` with the `'esr'` option.

## Implementation structure

- Define the magnet, spins, Zeeman interactions, temperature, coordinates, relaxation, basis, and propagator tolerance.
- Create the Spinach spin system and set the proton detection state and XiX experiment parameters.
- Run the steady-state calculation, plot the real proton `Lz` expectation value against microwave resonance offset, and save `xix_w_field_profile_single.fig`.
