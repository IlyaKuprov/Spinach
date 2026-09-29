# examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1_r.m)

`tppm_q_rep_time_ensemble_b1_r()` takes no inputs. It asks how the steady-state proton longitudinal-polarisation expectation varies with TPPM DNP repetition time after averaging over both electron–proton distance and microwave B1. The source header estimates hours of calculation time; this is the most expensive of the four repetition-time variants because it combines a four-node radial quadrature with a six-node B1 quadrature.

## Model and sequence

The source labels the model Q-band TPPM DNP, with `sys.magnet=1.2142`, spin temperature 80, and spins `{'E','1H'}` (the source gives no units for the magnet or temperature fields). The trityl g principal values are `[2.00319 2.00319 2.00258]`; the proton Zeeman entry is `[0 0 5]`, described in the source as a ppm guess. Euler-angle values `[0 10 0]` and `[0 0 10]` are converted from degrees to radians. The basis is `sphten-liouv` with no approximation, propagator chopping `1e-12`, and hygiene disabled.

The source obtains four distance nodes and weights from `gaussleg(3.5,20,3)` in Å and six B1 nodes and weights from `gaussleg(25e6,35e6,5)` in Hz. At each distance, the coordinates put the proton at `[0 0 r(n)]` Å relative to the electron. Relaxation uses `t1_t2`; its distance- and orientation-dependent nuclear R1 handle calls `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`. The source sets `inter.r1_rates={1e3 r1n_rate}`, `inter.r2_rates={200e3 50e3}`, diagonal relaxation retention, and `dibari` equilibrium; it does not annotate rate units.

TPPM settings are spins `{'E','1H'}`, grid `rep_2ang_800pts_sph`, 16 ns pulse duration, 300 loops, 120° second-pulse phase, −13 MHz added shift, and +2 MHz electron offset. The repetition-time scan is 30 log-spaced points from `10^-5` to `10^-2.7` s; the shot spacing is `rep_time - 2*nloops*pulse_dur`.

## Calculation and output

The detected operator is proton `Lz`, set with `state(spin_system,'Lz','1H')`. Each distance/B1/repetition-time point is evaluated by `powder(spin_system,@xixdnp_steady,localpar,'esr')`. Results are averaged over B1 with its quadrature weights and over distance with the source's explicit `r^2,wr` weighting (radial Jacobian), normalised by `sum(r.^2.*wr)`. The real proton expectation is plotted against repetition time in ms and saved as `tppm_q_rep_time_ensemble_b1_r.fig` in the current working directory. There is no declared function output.
