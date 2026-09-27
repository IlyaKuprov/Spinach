# examples/dnp_sol/steady_state/xix_q_rep_time_single.m

- Signature: `xix_q_rep_time_single()`

## Purpose

Simulate a steady-state XiX DNP repetition-time scan at Q-band; the source estimates a calculation time of seconds.

## Physical / mathematical content

- An electron–proton system at 80 K has specified Zeeman interactions and coordinates separated by 3.5 Å. The simulation includes electron and proton relaxation, an inverted-phase second XiX pulse, and powder averaging; it detects the proton $L_z$ expectation value.

## Numerical / algorithmic content

- The code scans 30 logarithmically spaced repetition times from $10^{-5}$ to $10^{-3}$ seconds. For each time, it subtracts the duration of 36 two-pulse XiX blocks to obtain the shot spacing and runs `powder(spin_system,@xixdnp_steady,localpar,'esr')` in a `parfor` loop. It plots the real signal against repetition time in milliseconds and saves `xix_q_rep_time_single.fig`.

## Implementation structure

- Set the Q-band field, electron and proton Zeeman interactions, temperature, and Cartesian coordinates; obtain the electron–nuclear distance.
- Specify the unrestricted spherical-tensor Liouville-space basis, propagator tolerance, and hygiene option.
- Configure distance- and orientation-dependent proton longitudinal relaxation, other relaxation rates, and equilibrium; create the spin system and proton detection operator.
- Set the powder grid, irradiation power, pulse duration, XiX block count, pulse phase, and frequency offsets before scanning repetition times.
