# examples/dnp_sol/steady_state/top_q_rep_time_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_rep_time_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_rep_time_single.m)

- Signature: `top_q_rep_time_single()`; source runtime estimate: seconds.

## Question and choice

How does TOP DNP steady-state proton longitudinal polarisation vary with repetition time for one electron–proton pair and one fixed microwave nutation frequency? This baseline uses 3.5 Å distance and 18 MHz B1, with no ensemble averaging.

## Model and scan

No arguments. Source settings: `sys.magnet=1.2142` (labelled Q-band), E/1H isotopes, trityl g values `[2.00319 2.00319 2.00258]`, proton Zeeman values `[0 0 5]` (source comment: ppm guess), Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, spin temperature 80, and coordinates (0,0,0)/(0,0,3.5); source derives `r_en` from the z-coordinate. Full `sphten-liouv` basis/no approximation; chopping tolerance `1e-12`; hygiene disabled. B1 is `parameters.irr_powers=18e6` Hz. Repetition time scans 30 log-spaced points from `1e-5` to `1e-3` seconds (10 μs–1 ms).

TOP settings: E/1H; spherical grid `rep_2ang_800pts_sph`; 10 ns pulse; 14 ns delay; 300 blocks; `addshift=-13e6`; `el_offs=95e6`. Shot spacing is repetition time minus `300*(10 ns+14 ns)`, or 7.2 μs.

## Relaxation and calculation

Uses `inter.relaxation={'t1_t2'}`; `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`; `inter.r1_rates={1e3 r1n_rate}`; `inter.r2_rates={200e3 50e3}`; diagonal retention; `dibari` equilibrium. Rate units are not annotated. Proton detection is `coil_state(spin_system,'Lz','1H','exact')`; each scan point calls `powder(spin_system,@topdnp_steady,localpar,'esr')`.

## Output and limits

Plots real proton `I_z` expectation against repetition time in ms and saves `top_q_rep_time_single.fig` in the MATLAB current folder, not a separate numerical data file. TOP steady-state dynamics are delegated to `topdnp_steady`.
