# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1.m)

## Purpose

The source header estimates the calculation time as hours.

`xix_w_pulse_dur_ensemble_b1()` maps steady-state XiX DNP proton `I_z` over microwave offset and electron pulse duration, averaging over an electron Rabi-frequency (B1) quadrature. Distance is fixed at the single coordinate value `3.5`; there is no distance ensemble.

## Model and sequence

The model uses the same W-band electron–proton pair and spin parameters as the companion single-profile example: magnet setting `3.4`, isotopes `E` and `1H`, trityl g principal values `[2.00319 2.00319 2.00258]`, proton shift guess `[0 0 5]` ppm, Euler angles `[0 10 0]` and `[0 0 10]` degrees, and spin temperature `80`. It uses the full `sphten-liouv` basis, `prop_chop=1e-12`, `t1_t2` relaxation with `r1n_dnp`, R1 entries `1e3`, R2 entries `200e3` and `50e3`, diagonal retention, and Di Bari equilibrium.

The sequence sets the second pulse phase to `pi`, the additive shift to `-33e6`, and the repetition spacing per scan point to `167e-6 - 360e-9`. For each pulse duration, the XiX block count is `round(360e-9/(2*pulse_dur))`, using the source's 360 ns contact-time numerator. The code does not annotate units for the magnet setting, temperature, coordinate, relaxation-rate entries, shot-spacing field, or additive shift; the values are therefore shown as assigned rather than converted.

## Numerical scan and output

The pulse-duration vector contains 200 values from `2e-9` to `21e-9 s`. The offset vector contains 101 values from `-230e6` to `205e6 Hz` (the figure axis is MHz). The six-node B1 quadrature is defined by `[b1,wb1]=gaussleg(10e6,20e6,5)`; its inputs are labelled Hz in the source. For each B1 node, the script sets `irr_powers=b1(k)` and uses `parfor` over durations to call `powder(...,@xixdnp_steady,...,'esr')` on `rep_2ang_800pts_sph`. It then combines the B1-node results using `wb1` and plots the real proton expectation as an offset-by-duration map, with offset in MHz and duration in ns. The figure is saved as `xix_w_pulse_dur_ensemble_b1.fig`.

## Dependencies and scope

The script depends on Spinach system/basis/state/powder and plotting functions, `gaussleg`, `r1n_dnp`, and `xixdnp_steady`. The ensemble average is over B1 only; the distance remains fixed.
