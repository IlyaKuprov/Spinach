# examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_r.m)

## Purpose

Run the no-argument MATLAB function `xix_w_field_profile_ensemble_r()` with Spinach and its example helpers available; it calculates a steady-state XiX DNP microwave field profile averaged over an electron–proton distance distribution. The source labels the magnet W-band and sets `sys.magnet=3.4`; its calculation-time estimate is minutes.

## Model and scan

The model contains E and 1H with a trityl g-tensor [2.00319 2.00319 2.00258], proton shift [0 0 5] (the source labels the proton value a ppm guess), Euler angles [0 10 0] and [0 0 10] degrees, and spin-temperature value 80. A four-point Gauss–Legendre distance quadrature uses bounds 3.5 and 20; unlike the companion B1-and-distance source, this file does not annotate the distance unit. At each distance the coordinates are updated and `r1n_dnp` provides orientation- and distance-dependent proton longitudinal relaxation using additional parameters `2.00230`, `1e-3`, and `52`; the source sets `r1_rates={1e3 r1n_rate}`, `r2_rates={200e3 50e3}`, `t1_t2`, `rlx_keep='diagonal'`, and `equilibrium='dibari'`. The basis is `sphten-liouv` without approximation, with propagator chop tolerance `1e-12`.

The electron nutation frequency is fixed at 20e6 Hz; this variant averages distance only, not a B1 ensemble. At every distance it evaluates a 201-point offset scan from −300e6 to 300e6 Hz using `powder(spin_system,@xixdnp_steady,parameters,'esr')`. The XiX train uses 18 ns pulses, 10 blocks, an inverted second-pulse phase ( `pi` ), shot spacing 167e−6, additional shift −33e6, and powder grid `rep_2ang_800pts_sph`. The distance profiles are combined with Gauss–Legendre weights multiplied by `r^2`, the radial Jacobian, then normalised by the weighted sum.

## Output and dependencies

The real part of the distance-averaged DNP value is plotted as the proton `Lz` expectation value against microwave offset in MHz and saved to `xix_w_field_profile_ensemble_r.fig`. The driver uses Spinach's system, basis, detection, and powder functions and calls `gaussleg`, `r1n_dnp`, and `xixdnp_steady`. It saves the plot as a MATLAB figure, not a separate numeric result file.
