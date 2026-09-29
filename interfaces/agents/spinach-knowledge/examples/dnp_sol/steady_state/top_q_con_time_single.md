# examples/dnp_sol/steady_state/top_q_con_time_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_con_time_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_con_time_single.m)

- Signature: `top_q_con_time_single()`

## Question and model

How does the steady-state proton longitudinal polarisation vary with TOP contact time for two irradiation conditions at one fixed electron–proton separation? This is the fixed-distance counterpart to the distance-ensemble contact-time scripts: the coordinates place the electron and proton 3.5 Å apart. The model uses a Q-band magnet setting of 1.2142, trityl electron g principal values `[2.00319 2.00319 2.00258]`, proton shift values `[0 0 5]`, Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, and spin temperature 80 K.

## Scan and sequence settings

The scan is `nloops=1:256` TOP blocks, each with a 10 ns pulse and 14 ns delay, so the plotted contact time is 24 ns times the loop count. Setting A uses 18 MHz irradiation and 95 MHz electron offset; setting B uses 33 MHz and 92 MHz. Their shot spacings are 102 μs and 153 μs, respectively, less the full pulse-train duration. The proton detector is `state(spin_system,'Lz','1H')`. The experiment uses spins `E` and `1H`, grid `rep_2ang_800pts_sph`, and `addshift=-13e6`; the source assigns R1 entries `1e3` and R2 values `200e3` and `50e3` (units are not annotated). `hygiene` is disabled.

For each loop count and setting, the steady state is obtained with `powder(spin_system,@topdnp_steady,localpar,'esr')`. The orientation-dependent proton R1 function `r1n_dnp` uses the fixed distance 3.5 Å. The relaxation model is `t1_t2`, with diagonal terms retained and `dibari` equilibrium; the basis is `sphten-liouv` with no approximation and propagator chopping tolerance `1e-12`.

## Output and limits

The figure compares the real proton `I_z` expectation for the two settings against total contact time and is saved as `top_q_con_time_single.fig`. There is no distance or B1 averaging in this script. The source estimates hours of calculation and saves a figure rather than a numeric results table.
