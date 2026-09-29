# examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_r.m)

`tppm_q_rep_time_ensemble_r()` takes no inputs. It asks how the steady-state proton longitudinal-polarisation expectation varies with TPPM DNP repetition time when electron–proton distance is averaged, at one fixed microwave nutation frequency. Unlike the `ensemble_b1` page, this source has no B1 quadrature: it fixes `parameters.irr_powers=33e6` Hz. The source header estimates minutes of calculation time.

## Model and sequence

The source labels the simulation Q-band TPPM DNP and sets `sys.magnet=1.2142`, spin temperature 80, and spins `{'E','1H'}`; units are not annotated for the magnet and temperature values. The trityl g principal values are `[2.00319 2.00319 2.00258]`; the proton Zeeman entry `[0 0 5]` is commented as a ppm guess. Euler angles `[0 10 0]` and `[0 0 10]` are converted from degrees to radians. It uses a `sphten-liouv` basis without approximation, propagator chopping `1e-12`, and disables hygiene.

Distance is represented by four Gauss–Legendre nodes and weights from `gaussleg(3.5,20,3)` (Å); for each node, the proton coordinate is `[0 0 r(n)]` Å. Relaxation uses `t1_t2`; the nuclear R1 handle calls `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`. The other source settings are `inter.r1_rates={1e3 r1n_rate}`, `inter.r2_rates={200e3 50e3}`, diagonal retention, and `dibari` equilibrium (rate units are not stated).

The TPPM parameters are spins `{'E','1H'}`, orientation grid `rep_2ang_800pts_sph`, fixed electron nutation frequency 33 MHz, pulse duration 16 ns, 300 loops, second-pulse phase 120°, added shift −13 MHz, and electron offset +2 MHz. It scans 30 logarithmically spaced repetition times from `10^-5` to `10^-2.7` s and sets shot spacing to `rep_time - 2*nloops*pulse_dur`.

## Calculation and output

The detection operator is proton `Lz`, created by `state(spin_system,'Lz','1H')`. For each distance and repetition time, the source calls `powder(spin_system,@xixdnp_steady,localpar,'esr')`. It combines distance points using the explicitly coded `r^2 * wr` radial weighting and normalises by `sum(r.^2.*wr)`. The real proton expectation is plotted against repetition time in ms and saved as `tppm_q_rep_time_ensemble_r.fig` in the current working directory. The function has no declared output argument.
