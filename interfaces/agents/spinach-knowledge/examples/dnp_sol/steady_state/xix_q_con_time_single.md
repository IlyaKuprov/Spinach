# examples/dnp_sol/steady_state/xix_q_con_time_single.m

Signature: `xix_q_con_time_single()`

This is the single-geometry XiX DNP contact-time calculation at steady state. Unlike the three ensemble variants, it fixes the electron–proton separation at 3.5 Å and does not average over distance. The source estimates the calculation time as minutes.

## Setup and scan

The E–¹H model uses `sys.magnet=1.2142` (Q-band setting), spin temperature `80`, trityl g values `[2.00319 2.00319 2.00258]`, proton Zeeman entry `[0 0 5]` (source comment: ppm guess), and Euler-angle entries `[0 10 0]` and `[0 0 10]` degrees. The coordinates place the proton on the z axis at 3.500; that distance is passed to `r1n_dnp`, whose source documents its `r` argument in Angstrom. The source derives the scalar electron–nuclear separation from those coordinates and calls `r1n_dnp` for the orientation-dependent nuclear R1, using arguments `sys.magnet`, temperature, `2.00230`, `1e-3`, `52`, that separation, and `bet`. It sets `inter.r1_rates={1e3,r1n_rate}` and `inter.r2_rates={200e3,50e3}`, with `t1_t2` relaxation, diagonal retention, and Di Bari equilibrium.

The full `sphten-liouv` basis has no approximation; propagator chopping tolerance is `1e-12`. The contact scan is XiX loop counts 1–64, with two 48 ns pulses per loop. The source uses `phase=pi` (inverted second pulse), `18e6` electron nutation frequency (Hz), grid `rep_2ang_800pts_sph`, `addshift=-13e6`, and `el_offs=61e6`; shot spacing is 153 μs minus the pulse-train duration. Contact time is twice pulse duration times loop count and is plotted in μs.

## Run dependencies and output

Requires Spinach MATLAB functions including `powder` and system/basis/state and plotting routines, plus [`r1n_dnp`](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r1n_dnp.m) and [`xixdnp_steady`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m). The mapped source is [`xix_q_con_time_single.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_single.m).

The no-argument function returns no MATLAB output. It plots the real proton `Lz` expectation value against contact time and saves `xix_q_con_time_single.fig` in the current directory. This is one fixed-separation result, not a distance-distribution average; the source does not emit a table or array of values as a function output.
