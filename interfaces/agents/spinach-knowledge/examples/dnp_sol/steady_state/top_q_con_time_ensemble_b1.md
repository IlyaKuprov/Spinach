# examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1.m

## Purpose and call

Call `top_q_con_time_ensemble_b1()` with no arguments from MATLAB with Spinach and the example helpers on the path. The function compares steady-state TOP DNP proton polarisation versus total contact time for two separate electron-Rabi-frequency ensembles. It declares no return value and saves a figure; the source estimates hours of calculation.

Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1.m).

## Spin model and relaxation

The system is an `E`/`1H` pair at a Q-band magnet setting `sys.magnet=1.2142`. The trityl electron Zeeman principal values are `[2.00319 2.00319 2.00258]`, proton values `[0 0 5]` (called a ppm guess), and orientations `(pi/180)*{[0 10 0],[0 0 10]}`. The spin-temperature parameter is `80`; no unit is stated. The pair is separated by `3.500` along z. The basis is `sphten-liouv`, with no approximation, `prop_chop=1e-12`, and hygiene disabled.

Relaxation is `t1_t2`, diagonal retention, and DiBari equilibrium. The orientation- and distance-dependent nuclear R1 function is called as `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`. Other configured relaxation entries are `1e3` and `[200e3 50e3]`; the source does not annotate units for these values. See [r1n_dnp.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r1n_dnp.m).

## TOP contact-time and B1 ensembles

The experiment uses `pulse_dur=10e-9` seconds, `delay_dur=14e-9`, `addshift=-13e6`, and powder grid `rep_2ang_800pts_sph`. Detection uses `state(spin_system,'Lz','1H')` and the experiment spins are `{'E','1H'}`. It scans `loop_counts=1:256`; contact time is calculated as `(pulse_dur+delay_dur)*loop_counts` and plotted in microseconds. For each loop count and B1 node, `powder` calls `topdnp_steady` in the `esr` context. Each ensemble has six Gauss–Legendre nodes: 10–20 MHz with `el_offs=95e6` (A), and 25–35 MHz with `el_offs=92e6` (B). The source legend names the curves “TOP, 15 MHz” and “TOP, 30 MHz”; those labels accompany the respective finite B1 ranges. Each B1 average uses quadrature weights normalised by their sum. See [topdnp_steady.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/topdnp_steady.m).

Shot spacing is set separately as `102e-6 - pulses_dur` for A and `153e-6 - pulses_dur` for B. These are the literal source expressions; units are not annotated in those assignments.

## Dependencies

The source uses Spinach system construction and propagation functions `create`, `basis`, `state`, and `powder`, the `gaussleg` quadrature helper, and MATLAB `parfor`; plotting uses `kfigure`, `kylabel`, and `klegend` before MATLAB `savefig`. The TOP propagator and R1 helper are linked above.

## Output and scope

The saved figure `top_q_con_time_ensemble_b1.fig` plots the real proton longitudinal expectation value against total contact time for the two B1 ensembles. It represents the specified pair, TOP helper, relaxation model, and six-node quadratures; the source provides no numerical result array as a function return.
