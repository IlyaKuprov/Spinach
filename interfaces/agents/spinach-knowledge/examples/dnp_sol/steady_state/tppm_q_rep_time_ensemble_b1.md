# examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_b1.m)

`tppm_q_rep_time_ensemble_b1()` takes no inputs. It asks how the steady-state proton longitudinal-polarisation expectation changes with DNP repetition time when microwave B1 is averaged over a six-node distribution, while the electron–proton separation is fixed at 3.5 Å. The source header estimates minutes of calculation time.

## Model and sequence

The source labels this a Q-band TPPM DNP simulation. It sets `sys.magnet=1.2142` and a spin temperature of 80 (the source does not annotate units for either field). The spins are `{'E','1H'}`; the trityl g principal values are `[2.00319 2.00319 2.00258]`, and the proton Zeeman entry is `[0 0 5]` (commented as a ppm guess). Euler angles are supplied as degree values converted to radians: `[0 10 0]` and `[0 0 10]`. The coordinates place the proton 3.500 Å along z from the electron.

The basis is `sphten-liouv` with `approximation='none'`; propagator chopping is `1e-12`. Relaxation uses `t1_t2`, diagonal retention and `dibari` equilibrium. The nuclear R1 function handle calls `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`; the other listed rates are `inter.r1_rates={1e3 r1n_rate}` and `inter.r2_rates={200e3 50e3}` (units are not stated in this example).

For the TPPM train, the source sets spins `{'E','1H'}`, orientation grid `rep_2ang_800pts_sph`, 16 ns pulse duration, 300 loops, second-pulse phase 120°, added shift −13 MHz, and electron offset +2 MHz. The six B1 nodes and weights come from `gaussleg(25e6,35e6,5)` (Hz). The repetition-time axis has 30 logarithmic points from `10^-5` to `10^-2.7` s. For each point, shot spacing is set to `rep_time - 2*nloops*pulse_dur`.

## Calculation and output

The proton detection operator is `state(spin_system,'Lz','1H')`. For each B1 node and repetition time, the source calls `powder(spin_system,@xixdnp_steady,localpar,'esr')`; it then averages the result with the B1 quadrature weights. The plot uses the real part of the proton `Lz` expectation against repetition time in ms, and is saved in the current working directory as `tppm_q_rep_time_ensemble_b1.fig`. The function has no declared output argument; the plotted data are local to the function.
