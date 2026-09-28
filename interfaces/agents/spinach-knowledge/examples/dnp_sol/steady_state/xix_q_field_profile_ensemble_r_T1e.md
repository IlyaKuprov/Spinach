# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1e.m

- Signature: `xix_q_field_profile_ensemble_r_T1e()`

## Purpose

Simulates the dependence of steady-state XiX DNP field profiles on electron longitudinal relaxation time, averaging over an electron–proton distance ensemble. The source estimates a calculation time of minutes.

## Physical / mathematical content

- Models an electron–proton pair at a 1.2142 T Q-band magnetic field and 80 K, with a trityl electron g-tensor and a proton Zeeman-shift estimate.
- Sweeps electron T1 values of 10, 3.0, 1.0, 0.3, and 0.1 ms. For each distance, the proton longitudinal relaxation rate is supplied by an orientation-dependent `r1n_dnp` function handle; electron R1 is `1/T1e`.
- Averages the calculated proton steady-state signal over electron–proton distances using Gauss–Legendre weights and an `r²` radial Jacobian, normalized by the sum of the weighted `r²` values.

## Numerical / algorithmic content

- Uses three Gauss–Legendre distance points from 3.5 to 20 and 201 microwave offsets from −100 to 100 MHz.
- For each distance, creates a Spinach electron–proton system in an unrestricted `sphten-liouv` basis, detects proton `Lz`, and evaluates `powder(spin_system,@xixdnp_steady,parameters,'esr')` on the `rep_2ang_800pts_sph` grid.
- Sets an 18 MHz electron nutation frequency, 48 ns pulse duration, 36 XiX blocks, an inverted second-pulse phase, a −13 MHz additional shift, and shot spacing calculated from 204 µs minus the two pulses per block.

## Implementation structure

- Initializes a figure, runs the field-profile calculation for each T1e value, and plots the real, distance-averaged proton signal against microwave offset in MHz.
- Adds a legend for the five T1e values and saves the figure as `xix_q_field_profile_ensemble_r_T1e.fig`.