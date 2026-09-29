# examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1_r.m)

- Signature: `top_q_rep_time_ensemble_b1_r()`; source runtime estimate: hours.

## Question and choice

How does TOP DNP steady-state proton longitudinal polarisation vary with repetition time after averaging over both electron–proton distance and microwave B1? This combined ensemble differs from the B1-only and distance-only variants by varying both axes.

## Model and scan

No arguments. Source settings: `sys.magnet=1.2142` (labelled Q-band), E/1H isotopes, trityl g values `[2.00319 2.00319 2.00258]`, proton Zeeman values `[0 0 5]` (source comment: ppm guess), Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, spin temperature 80, full `sphten-liouv` basis/no approximation, chopping tolerance `1e-12`, hygiene disabled. Distance quadrature is `[r,wr]=gaussleg(3.5,20,3)` (source comment: Å); B1 quadrature is `[b1,wb1]=gaussleg(10e6,20e6,5)` (Hz). Source does not further specify the quadrature argument convention. Coordinates vary as (0,0,0)/(0,0,r(n)). Repetition time scans 30 log-spaced points from `1e-5` to `1e-3` seconds (10 μs–1 ms).

For each B1 node, electron nutation frequency is `b1(k)`. TOP settings: E/1H; spherical grid `rep_2ang_800pts_sph`; 10 ns pulse; 14 ns delay; 300 blocks; `addshift=-13e6`; `el_offs=95e6`. Shot spacing subtracts the 7.2 μs pulse-train duration from repetition time.

## Relaxation, averaging, and output

At each distance, source uses `inter.relaxation={'t1_t2'}` and `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`, rates `inter.r1_rates={1e3 r1n_rate}` and `inter.r2_rates={200e3 50e3}`, diagonal retention, and `dibari` equilibrium (rate units not annotated). Proton `Lz` is detected. Each point calls `powder(spin_system,@topdnp_steady,localpar,'esr')`. It normalises the B1-weighted average by `sum(wb1)`, then distance-averages using `r^2*wr` (explicitly identified as the Jacobian), normalised by `sum(r.^2.*wr)`. The real proton `I_z` expectation is plotted against repetition time in ms; saves `top_q_rep_time_ensemble_b1_r.fig` in the MATLAB current folder, not separate data. TOP steady-state dynamics are delegated to `topdnp_steady`.
