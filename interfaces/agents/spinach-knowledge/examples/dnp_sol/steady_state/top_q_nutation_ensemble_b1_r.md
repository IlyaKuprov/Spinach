# examples/dnp_sol/steady_state/top_q_nutation_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_nutation_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_nutation_ensemble_b1_r.m)

- Signature: `top_q_nutation_ensemble_b1_r()`

## Question and model

How does the steady-state TOP DNP proton-polarisation field profile over microwave offset change across selected electron nutation frequencies, when distance and B1 are distributed? The outer script evaluates `nu= [6.8 9.6 13.5 17.5 25 36]` MHz. Each profile uses the Q-band magnet setting 1.2142, trityl electron g principal values `[2.00319 2.00319 2.00258]`, proton shift values `[0 0 5]`, Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, and spin temperature 80 K.

## Offset sweep and averaging

For each `nu`, the local helper `top_field_profile_b1_r(nu)` scans 10 microwave resonance offsets linearly from 88 to 97 MHz. It calls `gaussleg` with orders 3 and 5, which return four distance nodes from 3.5 to 20 Å and six B1 nodes from `0.2*nu` to `1.2*nu` Hz: 24 weighted distance/B1 combinations per offset. The TOP sequence has 300 blocks, each with a 10 ns pulse and 14 ns delay; shot spacing is 153 μs minus the pulse-train duration. The proton detector is `coil_state(spin_system,'Lz','1H','exact')`. Orientation-dependent proton R1 is set through `r1n_dnp`. The source assigns R1 entries `1e3` and R2 values `200e3` and `50e3` (units are not annotated), uses spins `E` and `1H` with grid `rep_2ang_800pts_sph`, sets `addshift=-13e6`, and disables `hygiene`.

Every offset/B1 point calls `powder(spin_system,@topdnp_steady,parameters,'esr')`. The relaxation model is `t1_t2`, with diagonal terms retained and `dibari` equilibrium; the basis is `sphten-liouv` without approximation and propagator chopping tolerance `1e-12`.

## Output and limits

The result is averaged over B1 weights and then distance weights with radial Jacobian `r^2`. Each curve plots the real proton `I_z` expectation against microwave offset; the 3-D axes are nutation frequency, offset, and steady-state expectation. The figure is saved as `top_q_nutation_ensemble_b1_r.fig`. The source estimates minutes of calculation; no numeric results table is written.
