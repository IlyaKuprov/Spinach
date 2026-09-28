# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2e.m

- Signature: `xix_q_field_profile_ensemble_r_T2e()`

## Purpose

Simulate steady-state XiX DNP field profiles for five electron transverse relaxation times, averaging each profile over an electron–proton distance ensemble. The source estimates a calculation time of minutes.

## Physical / mathematical content

- Models an electron and a proton at a Q-band magnetic field of 1.2142 T and a spin temperature of 80 K, using a trityl electron g-tensor and an estimated proton chemical shift.
- Evaluates $T_{2e}$ values of 50, 15, 5, 1.5, and 0.5 μs. The electron transverse relaxation rate is set to $1/T_{2e}$; the proton transverse rate is 50 kHz. The proton longitudinal relaxation rate depends on distance and orientation through `r1n_dnp`.
- Calculates the steady-state proton $L_z$ response across electron microwave resonance offsets from −100 to 100 MHz. Distance-ensemble averaging uses Gauss–Legendre weights and an $r^2$ Jacobian, normalized by the total weighted $r^2$.

## Numerical / algorithmic content

- Uses three Gauss–Legendre distance points from 3.5 to 20 and 201 microwave offsets. For each distance, it builds a Spinach spin system in an untruncated spherical-tensor Liouville basis and calls `powder(spin_system,@xixdnp_steady,parameters,'esr')` on the `rep_2ang_800pts_sph` orientation grid.
- Sets an electron nutation frequency of 18 MHz, 48 ns pulse duration, 36 XiX blocks, an inverted second-pulse phase, and a −13 MHz additional shift. Shot spacing is calculated as `204e-6 - 2*parameters.nloops*parameters.pulse_dur`.

## Implementation structure

- The main function initializes the figure, calls the profile helper once per $T_{2e}$ value, adds a legend, and saves `xix_q_field_profile_ensemble_r_T2e.fig`.
- The helper computes the distance-resolved steady-state profiles, averages them over distance, and plots the real response against microwave offset in MHz.