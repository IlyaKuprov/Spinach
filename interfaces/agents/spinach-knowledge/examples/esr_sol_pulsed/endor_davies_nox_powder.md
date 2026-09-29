# examples/esr_sol_pulsed/endor_davies_nox_powder.m

Call `endor_davies_nox_powder()` with no arguments. This powder Davies ENDOR example uses an electron–14N nitroxide at 3.5 T; the source describes explicit soft pulses and orientation selection using a large spherical grid, evaluated with Fokker–Planck formalism. It is a brute-force time-domain calculation, with source-estimated runtime of hours.

## Spin system and simulation

- The electron g matrix is diagonal: 2.01045, 2.00641, 2.00211. The electron–14N coupling matrix is `[1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230] × 10^7` as specified in the source. The matrix units are not labelled there.
- The basis is `sphten-liouv` with no approximation; relaxation uses `t1_t2`, diagonal retention and zero equilibrium. The source supplies `r1_rates` and `r2_rates` as `{20e6, 0.5e6}`. Trajectory-level SSR is disabled.
- `powder` calls `@endor_davies` in the `esr` context on `rep_2ang_12800pts_sph`. The electron pulse is rank 2, phase π/2, frequency −300 MHz and duration 10 ns; the nuclear pulse is rank 3, phase π/2 and duration 100 ns. Nuclear frequency is swept from −80 to +80 MHz at 200 points. The source sets `e_pwr = 2π × 16.5 × 10^7` and `n_pwr = π × 10^7` (units are not annotated), and sets `offset = [−2×10^8, 0]`.
- The returned response is plotted as its real part against nuclear frequency in MHz; the vertical-axis label is `(RF on)/(RF off)`. The function displays a figure and does not save a data file.

## Use and limitations

This example is specifically the soft-pulse Davies variant: its orientation-selection effects and long run time distinguish it from the hard-pulse Mims nitroxide example. Its source does not provide a DOI.

Source: [examples/esr_sol_pulsed/endor_davies_nox_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/endor_davies_nox_powder.m).
