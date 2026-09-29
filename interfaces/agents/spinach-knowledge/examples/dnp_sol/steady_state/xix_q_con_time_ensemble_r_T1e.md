# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1e.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1e.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1e.m)

- Signature: `xix_q_con_time_ensemble_r_T1e()`

## Purpose

Compares steady-state XiX proton contact-time curves for five electron longitudinal-relaxation times while averaging over electron–proton distance. It fixes B1 at 18 MHz; it is not the B1-ensemble scan in `xix_q_con_time_ensemble_b1` or `xix_q_con_time_ensemble_b1_r`.

## Model and sweep

The outer function evaluates `T1e=[10e-3,3.0e-3,1.0e-3,0.3e-3,0.1e-3]` s, i.e. 10, 3.0, 1.0, 0.3, and 0.1 ms. For each value it calls the local helper `xix_contact_curve_ensemble_r(T1e)`. The helper builds the same Q-band electron–`1H` pair at `sys.magnet=1.2142` (1.2142 T), with electron Zeeman principal values `[2.00319 2.00319 2.00258]`, proton shift `[0 0 5]` (source: ppm guess), Euler-angle triplets `[0 10 0]` degrees converted to radians, and spin temperature 80 K. The four-node distance quadrature is `gaussleg(3.5,20,3)` Å, and coordinates are set to `[0 0 0]` and `[0 0 r]` at each node.

The basis is `sphten-liouv` with no approximation, propagator chop tolerance `1e-12`, and hygiene disabled. Relaxation uses `t1_t2`, diagonal terms, and `dibari` equilibrium. The electron R1 rate is `1/T1e`; electron R2 and proton R2 are 200000 and 50000. Proton R1 is supplied by `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`, a distance- and orientation-dependent function handle. Rate units are not annotated in the source.

## Contact-time calculation and output

For each T1e and distance node, the helper uses fixed 18 MHz electron nutation, grid `rep_2ang_800pts_sph`, 48 ns pulse duration, π second-pulse phase, `addshift=-13e6`, and `el_offs=61e6`. It scans `nloops=1:64` and sets shot spacing to 153 μs minus total pulse duration. Each point calls `powder(spin_system,@xixdnp_steady,localpar,'esr')`; the distance result is averaged using quadrature weights and the radial (r^2) Jacobian. Total contact time is `2*nloops*48e-9` s (96 ns to 6.144 μs).

The five real proton (L_z) expectation-value curves are overlaid and labelled by T1e. The figure limits are 0–6 μs horizontally and 0–1.7e-3 vertically; it is saved as `xix_q_con_time_ensemble_r_T1e.fig`. The function returns no explicit MATLAB output. The source comment estimates hours of calculation time, not a measured runtime here.
