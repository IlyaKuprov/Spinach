# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m

- Signature: `xix_q_field_profile_ensemble_r()`

## Purpose

Simulate a steady-state XiX DNP microwave-offset profile at Q-band, averaged over an electron–proton distance ensemble. The source estimates a calculation time of minutes.

## Physical / mathematical content

The system contains an electron and a proton at 80 K in a 1.2142 T field. It uses a trityl electron g-tensor, a proton Zeeman shift, and an electron–proton separation varied across the ensemble. The `t1_t2` relaxation model includes a proton longitudinal relaxation rate calculated by `r1n_dnp` as a function of distance and orientation. The final distance average uses the Gauss–Legendre weights and an `r²` Jacobian, normalized by the sum of those weights.

## Numerical / algorithmic content

Three Gauss–Legendre points span distances from 3.5 to 20. For each distance, the function creates a Spinach system in the `sphten-liouv` basis without basis approximation, detects proton `Lz`, and calls `powder(spin_system,@xixdnp_steady,parameters,'esr')`. The microwave resonance offsets comprise 201 points from −100 to 100 MHz. Powder averaging uses the `rep_2ang_800pts_sph` grid.

## Implementation structure

The XiX settings specify an 18 MHz electron nutation frequency, 48 ns pulse duration, 36 blocks, an inverted second-pulse phase of `pi`, and a −13 MHz additional shift. Shot spacing is calculated as `204e-6 - 2*nloops*pulse_dur`. After distance averaging, the function plots the real proton `Lz` expectation value against microwave resonance offset in MHz and saves `xix_q_field_profile_ensemble_r.fig`.