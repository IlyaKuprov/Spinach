# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2n.m

- Signature: `xix_q_field_profile_ensemble_r_T2n()`

## Purpose

Simulates how the nuclear transverse relaxation time T2n affects steady-state XiX DNP field profiles for an electron–proton distance ensemble. The source estimates a calculation time of minutes.

## Physical / mathematical content

- Models an electron and a proton at a Q-band magnetic field of 1.2142 T and a spin temperature of 80 K. The electron Zeeman tensor uses trityl g-values; the proton Zeeman tensor uses a ppm estimate.
- Sweeps T2n over 2000, 200, 20, 2, and 0.2 μs. For each value, the proton transverse relaxation rate is `1/T2n`; the proton longitudinal relaxation rate is calculated by `r1n_dnp` and depends on distance and orientation.
- Samples electron–proton separations from 3.5 to 20 using three Gauss–Legendre points. The resulting profiles are averaged with quadrature weights and an `r^2` Jacobian.

## Numerical / algorithmic content

- For each T2n and sampled distance, creates an electron–proton spin system in the unrestricted spherical-tensor Liouville-space basis, with `t1_t2` relaxation and a proton `Lz` detection state.
- Calculates steady-state XiX DNP with `powder(spin_system,@xixdnp_steady,parameters,'esr')` on the `rep_2ang_800pts_sph` grid. The microwave resonance offsets span −100 to 100 MHz in 201 points; the experiment uses an 18 MHz electron nutation frequency, 48 ns pulses, 36 XiX blocks, an inverted second-pulse phase, and a −13 MHz additional shift.
- Plots the real, distance-averaged proton steady-state signal against microwave resonance offset for each T2n, adds a legend, and saves `xix_q_field_profile_ensemble_r_T2n.fig`.
