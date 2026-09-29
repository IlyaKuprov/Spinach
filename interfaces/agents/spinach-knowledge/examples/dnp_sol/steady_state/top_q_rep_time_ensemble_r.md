# examples/dnp_sol/steady_state/top_q_rep_time_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_rep_time_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_rep_time_ensemble_r.m)

- Signature: `top_q_rep_time_ensemble_r()`; source runtime estimate: minutes.

## Question and choice

How does TOP DNP steady-state proton longitudinal polarisation vary with repetition time after averaging over electron–proton distance at fixed microwave nutation frequency? This is distance-only: unlike the B1 ensemble pages, it fixes B1 at 18 MHz.

## Model and scan

No arguments. Source settings: `sys.magnet=1.2142` (labelled Q-band), E/1H isotopes, trityl g values `[2.00319 2.00319 2.00258]`, proton Zeeman values `[0 0 5]` (source comment: ppm guess), Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, spin temperature 80, full `sphten-liouv` basis/no approximation, chopping tolerance `1e-12`, hygiene disabled. Distance quadrature is `[r,wr]=gaussleg(3.5,20,3)` (source comment: Å), with coordinates (0,0,0)/(0,0,r(n)); source does not further specify quadrature argument convention. B1 is fixed by `parameters.irr_powers=18e6` Hz. Repetition time is `logspace(-5,-3,30)` seconds (10 μs–1 ms).

TOP settings: E/1H; spherical grid `rep_2ang_800pts_sph`; 10 ns pulse; 14 ns delay; 300 blocks; `addshift=-13e6`; `el_offs=95e6`. Shot spacing subtracts the 7.2 μs pulse-train duration from each repetition time.

## Relaxation, averaging, and output

At each distance, uses `inter.relaxation={'t1_t2'}` and `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`; rates `inter.r1_rates={1e3 r1n_rate}` and `inter.r2_rates={200e3 50e3}`; diagonal retention and `dibari` equilibrium. Rate units are not annotated. Proton detection is `state(spin_system,'Lz','1H')`; each point calls `powder(spin_system,@topdnp_steady,localpar,'esr')`. Distance average uses weights and radial `r^2` Jacobian, normalised by `sum((r.^2).*wr)`. Plots real proton `I_z` expectation versus repetition time in ms and saves `top_q_rep_time_ensemble_r.fig` in the MATLAB current folder; no separate numerical data file. TOP steady-state dynamics are delegated to `topdnp_steady`.
