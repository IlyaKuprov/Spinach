# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1e.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1e.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1e.m)

- Signature: `xix_q_rep_time_ensemble_r_T1e()` (no input arguments).
- Source: [MATLAB implementation](../../../../../../examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1e.m).

## Purpose and variant

This variant scans electron longitudinal relaxation time while retaining the distance ensemble. It compares how the steady-state repetition-time profile changes over the source's five selected T1e values.

## Inputs, model, and scan

The function constructs an electron–proton system (`E`, `1H`) and plots the steady-state proton `I_Z` signal against XiX repetition time. The five values listed below define the overlaid traces; the internal repetition-time array is `logspace(-5,-3,30)` (30 points), plotted as `1e3*rep_time` on an axis labelled ms. For each trace it samples four electron–proton distances with `gaussleg(3.5,20,3)` (the source labels these limits in Å), calls `powder(spin_system,@xixdnp_steady,localpar,'esr')` at every distance and repetition time, then distance-averages using quadrature weights and the `r^2` radial Jacobian.

The shared spin/pulse setup is: `sys.magnet=1.2142`; Zeeman eigenvalues `[2.00319 2.00319 2.00258]` and `[0 0 5]`, Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`; spin temperature `80`; basis `sphten-liouv` with `approximation='none'`; and `sys.tols.prop_chop=1e-12`. The source does not state units for the magnet value or spin temperature. The relaxation model is `inter.relaxation={'t1_t2'}`, with `inter.rlx_keep='diagonal'` and `inter.equilibrium='dibari'`; hygiene is disabled with `sys.disable={'hygiene'}`. The XiX setup uses the `rep_2ang_800pts_sph` orientation grid, electron nutation frequency `18e6` Hz, pulse duration `48e-9` s, `36` XiX blocks, and `phase=pi` (the source comments that the second pulse has inverted phase). It sets `addshift=-13e6` and `el_offs=-39e6` without explicit units. Shot spacing is calculated as `rep_time - 2*nloops*pulse_dur`.

`T1e=[10e-3, 3.0e-3, 1.0e-3, 0.3e-3, 0.1e-3]` s. In the rate setup, electron R1 is `1/T1e`; the proton rate is supplied by `r1n_dnp(sys.magnet,inter.temperature,2.00230,T1e,52,r(n),bet)`, so that callback also receives the scanned T1e, distance, and powder polar angle. The fixed `inter.r2_rates` pair is `{200e3,50e3}`.

## Running, dependencies, and output

Run the zero-argument function with Spinach and the example helpers on the MATLAB path. It calls `create`, `basis`, `state`, `powder`, and `xixdnp_steady`; plotting uses `kfigure`, `kgrid`, `kxlabel`, `kylabel`, and `klegend`. The distance integration uses [`gaussleg.m`](../../../../../../kernel/grids/gaussleg.m), and the nuclear-relaxation callback uses [`r1n_dnp.m`](../../../../../../etc/textbook/r1n_dnp.m); no external parameter file is read. The function plots `-real(dnp)` and saves `xix_q_rep_time_ensemble_r_T1e.fig` in the current working directory; it does not write a separate numeric-results file. The source estimates calculation time in hours. The T1e sweep changes both electron R1 and the T1e argument used by the distance-/orientation-dependent proton-R1 callback; it is not only a change to electron relaxation.
