# examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_r.m

- Signature: `xix_w_field_profile_ensemble_r()`

## Purpose

Simulation of a steady-state XiX DNP field profile at a 3.4 T W-band magnet field, averaged over an electron–proton distance ensemble. Calculation time: minutes.

## Physical / mathematical content

- The model contains an electron and a proton with a trityl electron g-tensor, a proton Zeeman shift, and a temperature of 80 K. At each electron–proton distance, it uses `t1_t2` relaxation with an orientation- and distance-dependent proton longitudinal relaxation rate, and detects the proton `Lz` expectation value. The distance average includes an `r^2` Jacobian.

## Numerical / algorithmic content

- Three Gauss–Legendre points span electron–proton distances from 3.5 to 20. For each distance, the code constructs a spin system in an unrestricted `sphten-liouv` basis and calls `powder` with `@xixdnp_steady` over 201 microwave resonance offsets from −300 to 300 MHz. It combines the resulting profiles using the quadrature weights and `r^2` weighting.

## Implementation structure

- Set the magnet field, electron and proton isotopes, Zeeman interactions, temperature, basis, and propagator tolerance.
- Generate the distance quadrature and microwave-offset grid.
- For each distance, set Cartesian coordinates and relaxation, construct the spin system, and configure proton detection and XiX experiment parameters, including the microwave power, spherical grid, pulse duration, ten blocks, phase, shot spacing, and additional shift.
- Run the steady-state powder simulation, average over distance, plot the real proton expectation value against microwave offset, and save `xix_w_field_profile_ensemble_r.fig`.
