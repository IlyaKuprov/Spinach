# examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1.m)

- Signature: `top_q_rep_time_ensemble_b1()`; source runtime estimate: minutes.

## Question and choice

How does TOP DNP steady-state proton longitudinal polarisation vary with repetition time when microwave B1 is distributed but electron–proton distance is fixed? This variant uses a 3.5 Å separation and no distance ensemble.

## Model and scan

No arguments. Source settings: `sys.magnet=1.2142` (labelled Q-band), isotopes E/1H, trityl electron g values `[2.00319 2.00319 2.00258]`, proton Zeeman values `[0 0 5]` (source comment: ppm guess), Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, spin temperature 80, and coordinates (0,0,0)/(0,0,3.5) (distance in Å). Basis is full `sphten-liouv` with no approximation; propagator chopping tolerance `1e-12`; hygiene disabled.

B1 nodes/weights are `[b1,wb1]=gaussleg(10e6,20e6,5)` (Hz per source comment); distance is fixed. The 30-point repetition-time scan is `logspace(-5,-3,30)` seconds (10 μs–1 ms). It sets electron nutation frequency to `b1(k)` at each B1 node. TOP settings: spins E/1H; spherical grid `rep_2ang_800pts_sph`; 10 ns pulse; 14 ns delay; 300 blocks; `addshift=-13e6`; `el_offs=95e6`. Shot spacing subtracts 300 pulse-plus-delay periods (7.2 μs) from each repetition time.

## Relaxation and calculation

The source sets `inter.relaxation={'t1_t2'}`, callback `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`, `inter.r1_rates={1e3 r1n_rate}`, `inter.r2_rates={200e3 50e3}`, diagonal retention, and `dibari` equilibrium. Rate units are not annotated. Proton detection is `state(spin_system,'Lz','1H')`. Each point calls `powder(spin_system,@topdnp_steady,localpar,'esr')`; B1 results are averaged with normalised `wb1` weights.

## Output and limits

Plots the real part of proton `I_z` expectation against repetition time in ms and saves `top_q_rep_time_ensemble_b1.fig` in the MATLAB current folder; no separate numerical data file is saved. TOP steady-state dynamics are delegated to `topdnp_steady`.
