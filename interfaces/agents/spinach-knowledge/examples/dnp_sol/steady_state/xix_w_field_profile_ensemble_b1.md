# examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1.m)

## Purpose

Run the no-argument MATLAB function `xix_w_field_profile_ensemble_b1()` with Spinach and its example helpers available; it computes a steady-state XiX DNP field profile averaged over an electron nutation-frequency (B1) ensemble. The source labels the magnet setting W-band and assigns `sys.magnet=3.4`; its stated calculation-time estimate is minutes.

## Model and scan

The two-spin model uses isotopes E and 1H, the trityl g-tensor [2.00319 2.00319 2.00258], proton shift [0 0 5] (the source labels the proton value a ppm guess), Euler angles [0 10 0] and [0 0 10] degrees, and spin-temperature value 80. Coordinates place the electron and proton at z = 0 and 3.500; this source does not state a unit for that coordinate separation. Orientation- and distance-dependent proton longitudinal relaxation uses `r1n_dnp` with additional parameters `2.00230`, `1e-3`, and `52`; rates are set as `r1_rates={1e3 r1n_rate}` and `r2_rates={200e3 50e3}` within the `t1_t2` relaxation model; the basis uses `sphten-liouv` without approximation and a `1e-12` propagator chop tolerance.

A six-point Gauss–Legendre quadrature uses the B1 interval 10e6–20e6 Hz. At each point, the driver sets the electron nutation frequency and calls `powder(spin_system,@xixdnp_steady,parameters,'esr')` on 201 equally spaced microwave offsets from −300e6 to 300e6 Hz. The pulse train has 18 ns pulses, 10 XiX blocks, an inverted second-pulse phase ( `pi` ), shot spacing 167e−6, and additional shift −33e6; the orientation grid is `rep_2ang_800pts_sph`. The six offset profiles are combined using the Gauss–Legendre weights `wb1` and normalised by their sum. This is a B1 average only, not a distance ensemble.

## Output and dependencies

The real part of the averaged DNP values is plotted as the proton `Lz` expectation value against microwave offset in MHz and saved to `xix_w_field_profile_ensemble_b1.fig`. In addition to Spinach's system, basis, detection, and powder routines, the driver calls `gaussleg`, `r1n_dnp`, and `xixdnp_steady`. Its output is a figure; the source does not save a separate numeric result file.
