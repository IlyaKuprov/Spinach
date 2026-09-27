# examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1.m

- Signature: `xix_w_field_profile_ensemble_b1()`

## Purpose

Simulation of XiX DNP field profile in the steady state with electron Rabi frequency ensemble averaging. Calculation time: minutes.

## Physical / mathematical content

- Models an electron–proton pair at a 3.4 T W-band field and 80 K. The electron has a trityl g-tensor, the proton a specified chemical shift, and their Cartesian coordinates place them 3.5 distance units apart. Relaxation includes a distance- and orientation-dependent proton longitudinal rate supplied by `r1n_dnp`.
- Calculates the steady-state proton `Lz` expectation value across microwave offsets from −300 to 300 MHz for a ten-block XiX sequence with an inverted-phase second pulse.

## Numerical / algorithmic content

- Uses a full spherical-tensor Liouville-space basis, a spherical powder grid, and `xixdnp_steady` for the steady-state calculation. Five Gauss–Legendre points span electron nutation frequencies of 10–20 MHz; the resulting powder profiles are averaged with their quadrature weights.
- Plots the real, ensemble-averaged proton signal against microwave resonance offset and saves the figure as `xix_w_field_profile_ensemble_b1.fig`.

## Implementation structure

- Set the W-band magnet field and electron–proton isotopes.
- Specify Zeeman interactions, temperature, Cartesian coordinates, and electron–nuclear distance.
- Configure relaxation rates, equilibrium, the basis set, and propagator tolerance; create the Spinach spin system.
- Define proton detection, the B1 quadrature points, and XiX experiment parameters.
- Run the steady-state powder calculation for each B1 value, average the profiles, plot the result, and save the figure.
