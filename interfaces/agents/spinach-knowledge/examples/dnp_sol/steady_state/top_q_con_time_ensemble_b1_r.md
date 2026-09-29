# examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m)

- Signature: `top_q_con_time_ensemble_b1_r()`

## Question and model

How does the steady-state proton longitudinal polarisation change with TOP contact time for two irradiation settings when electron–proton distance and electron nutation frequency are both distributed? This Q-band electron–proton model uses `sys.magnet=1.2142`, trityl electron g principal values `[2.00319 2.00319 2.00258]`, proton shift values `[0 0 5]`, Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, and spin temperature 80 K.

## Scan and averaging

The contact-time axis is `nloops=1:256` TOP blocks, each comprising a 10 ns pulse and 14 ns delay (contact time = 24 ns × block count). Distance uses 3 Gauss–Legendre points over 3.5–20 Å. Two separate six-point B1 quadratures cover 10–20 MHz (A) and 25–35 MHz (B). For each distance, B1 point, and loop count, the code calls `powder(spin_system,@topdnp_steady,localpar,'esr')`. The proton coil is `state(spin_system,'Lz','1H')`; orientation-dependent proton R1 is supplied by `r1n_dnp` using the current distance and orientation angle `bet`. The source sets R1 entries to `1e3` and R2 values to `200e3` and `50e3` (units are not annotated), retains diagonal relaxation terms, and selects `dibari` equilibrium. The experiment uses spins `E` and `1H`, grid `rep_2ang_800pts_sph`, and `addshift=-13e6`.

Both settings use electron offset 95 MHz (A) or 92 MHz (B); their shot spacings are respectively 102 μs and 153 μs minus the pulse-train duration. The basis is `sphten-liouv` with no approximation, the propagator chopping tolerance is `1e-12`, and `hygiene` is disabled.

## Output and limits

B1 quadrature weights are applied first, then distance weights with the radial Jacobian `r^2`. The figure plots the real proton `I_z` expectation against total contact time for both ensembles and is saved as `top_q_con_time_ensemble_b1_r.fig`. The source estimates hours of calculation; it writes a figure, not a numeric results table.
