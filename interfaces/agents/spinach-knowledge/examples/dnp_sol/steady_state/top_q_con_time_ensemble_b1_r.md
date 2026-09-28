# examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1_r.m

- Signature: `top_q_con_time_ensemble_b1_r()`

## Purpose

Calculates the proton longitudinal expectation value versus total TOP DNP contact time, averaging over electron-proton distance and electron Rabi-frequency ensembles. The source comments estimate the calculation takes hours.

## Physical / mathematical content

The model is an electron-proton pair in a Q-band field (1.2142 T), with the same trityl electron g-tensor and 1H shift parameters as the companion fixed-distance example and spin temperature 80 K. The distance is sampled at three Gauss-Legendre nodes from 3.5 to 20 Å. For each sampled distance, the proton coordinate is set along z and the distance-dependent proton relaxation rate is evaluated with `r1n_dnp`. The full sphten-liouv basis, diagonal relaxation, dibari equilibrium, and `rep_2ang_800pts_sph` powder grid are used.

For each distance the script runs loop counts 1–256 for two five-node B1 distributions (10–20 MHz and 25–35 MHz), with 10 ns pulses and 14 ns delays. The two parameter sets use electron offsets of 95 MHz and 92 MHz and shot spacings of 102 μs and 153 μs, respectively, each reduced by the pulse-train duration. Each point is evaluated by `powder` with `@topdnp_steady` in `esr` mode. Results are averaged over B1 quadrature weights and then over distance with the radial Jacobian factor r². The plot compares proton Iz versus total contact time for the two ensembles and is saved as `top_q_con_time_ensemble_b1_r.fig`.
