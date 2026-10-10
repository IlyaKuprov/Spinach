# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2n.m

Signature: `xix_q_con_time_ensemble_r_T2n()`

This distance-ensemble XiX steady-state DNP example compares proton-contact curves while scanning nuclear T2. The source estimates the calculation time as hours.

## Setup and scan

The E–¹H pair uses `sys.magnet=1.2142` (Q-band setting), spin temperature `80`, trityl g values `[2.00319 2.00319 2.00258]`, and proton Zeeman entry `[0 0 5]` (source comment: ppm guess). Euler-angle entries are `[0 10 0]` and `[0 0 10]` degrees. It uses the full `sphten-liouv` basis without approximation, `prop_chop=1e-12`, `sys.disable={'hygiene'}`, diagonal relaxation retention, and Di Bari equilibrium.

The scan is nuclear T2 `[2000e-6 200e-6 20e-6 2e-6 0.2e-6]` seconds. At each of four Gauss–Legendre distances from 3.5–20 Å (`gaussleg(3.5,20,3)`), the z-axis geometry is evaluated with orientation-dependent nuclear R1 from `r1n_dnp`, called with the source arguments `sys.magnet`, temperature, `2.00230`, `1e-3`, `52`, the current radius, and `bet`. Relaxation entries are `inter.r1_rates={1e3,r1n_rate}` and `inter.r2_rates={200e3,1/T2n}`; the distance values are combined with quadrature weights and the radial `r^2` Jacobian.

Each curve scans XiX loop counts 1–64. Pulses are 48 ns, with `phase=pi` for the inverted second pulse; contact time is twice pulse duration times loop count and is plotted in μs. The source updates shot spacing as 153 μs minus the two-pulse train duration. Other fixed settings: electron nutation frequency `18e6` Hz, grid `rep_2ang_800pts_sph`, `addshift=-13e6`, and `el_offs=61e6` (the latter two are source numeric settings without a unit annotation there).

## Run dependencies and output

Requires Spinach MATLAB functions including `gaussleg`, `powder`, and system/basis/state and plotting routines; also requires [`r1n_dnp`](https://github.com/IlyaKuprov/Spinach/blob/main/etc/textbook/r1n_dnp.m) and the steady-state sequence [`xixdnp_steady`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m). The mapped source is [`xix_q_con_time_ensemble_r_T2n.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2n.m).

The no-argument function returns no MATLAB output. It plots the real proton `Lz` expectation value after distance averaging and saves `xix_q_con_time_ensemble_r_T2n.fig` in the current directory. Only nuclear T2 varies among the five curves; the electron R2 entry, nuclear R1 model, distance quadrature, and pulse protocol are shared.
