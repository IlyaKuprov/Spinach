# examples/dnp_sol/steady_state/tppm_q_rep_time_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_rep_time_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_rep_time_single.m)

`tppm_q_rep_time_single()` takes no inputs. It asks how the steady-state proton longitudinal-polarisation expectation changes with TPPM DNP repetition time for one fixed electron–proton separation and one fixed microwave nutation frequency. This is the single-condition baseline among the four variants: it performs no distance or B1 ensemble averaging. The source header estimates seconds of calculation time.

## Model and sequence

The source labels the model Q-band TPPM DNP. It sets `sys.magnet=1.2142` and spin temperature 80; no units are attached to those values in the source. The spins are `{'E','1H'}`; trityl g principal values are `[2.00319 2.00319 2.00258]`, and the proton Zeeman entry `[0 0 5]` is described as a ppm guess. Euler angles `[0 10 0]` and `[0 0 10]` are converted from degrees to radians. Coordinates fix the electron–proton separation at 3.5 Å along z. The basis is `sphten-liouv` with `approximation='none'`; propagator chopping is `1e-12`.

Relaxation is set to `t1_t2`, with the nuclear R1 function handle calling `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`. The source also sets `inter.r1_rates={1e3 r1n_rate}`, `inter.r2_rates={200e3 50e3}`, diagonal relaxation retention, and `dibari` equilibrium; rate units are not specified in the file.

The sequence parameters are spins `{'E','1H'}`, orientation grid `rep_2ang_800pts_sph`, fixed electron nutation frequency 33 MHz, 16 ns pulse duration, 300 loops, second-pulse phase 120°, added shift −13 MHz, and electron offset +2 MHz. The 30-point repetition-time axis is logarithmic from `10^-5` to `10^-2.7` s; each shot spacing is `rep_time - 2*nloops*pulse_dur`.

## Calculation and output

The proton detection operator is `coil_state(spin_system,'Lz','1H','exact')`. Each repetition-time point is evaluated with `powder(spin_system,@xixdnp_steady,localpar,'esr')`. The script plots the real proton `Lz` expectation value against repetition time in ms and saves `tppm_q_rep_time_single.fig` in the current working directory. The function declares no output argument.
