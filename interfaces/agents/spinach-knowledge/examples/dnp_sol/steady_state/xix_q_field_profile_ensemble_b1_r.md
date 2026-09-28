# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1_r.m

- Signature: `xix_q_field_profile_ensemble_b1_r()`

## Purpose

Simulate a steady-state XiX DNP field profile at Q band, averaging over electron–proton distance and electron Rabi frequency. The source estimates a calculation time of minutes.

## Physical / mathematical content

- Models an electron (`E`) and proton (`1H`) at a 1.2142 T magnetic field and 80 K, with a trityl electron g-tensor and an estimated proton chemical shift.
- Places the proton at each sampled distance from the electron. Uses `t1_t2` relaxation with a distance- and orientation-dependent proton longitudinal relaxation rate supplied by `r1n_dnp`; the equilibrium model is `dibari`.
- Detects proton `Lz` while sweeping microwave resonance offsets from −100 to 100 MHz.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism without basis approximation and a spherical powder grid (`rep_2ang_800pts_sph`). At each distance and electron nutation frequency, `powder` runs `xixdnp_steady` with the `esr` setting.
- Sets 48 ns pulses, 36 XiX DNP blocks, an inverted second-pulse phase (`pi`), a −13 MHz additional shift, and shot spacing calculated as `204e-6 - 2*nloops*pulse_dur` seconds.
- Samples distances from 3.5 to 20 Å at three Gauss–Legendre points and electron Rabi frequencies from 10 to 20 MHz at five points. Averages the resulting profiles using the B1 quadrature weights, then the distance quadrature weights multiplied by the radial factor `r^2`.

## Output

Plots the real proton `Lz` expectation value against microwave resonance offset in MHz and saves the figure as `xix_q_field_profile_ensemble_b1_r.fig`.