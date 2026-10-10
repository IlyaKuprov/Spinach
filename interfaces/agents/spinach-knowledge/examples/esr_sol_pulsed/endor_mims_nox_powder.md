# examples/esr_sol_pulsed/endor_mims_nox_powder.m

Call `endor_mims_nox_powder()` with no arguments. It simulates Mims ENDOR for a nitroxide radical powder with ideal hard pulses; the source estimates a runtime of seconds.

## Spin system and processing

- The system is electron plus 14N at 3.5 T. Its diagonal electron g values are 2.01045, 2.00641 and 2.00211. The electron–14N coupling matrix is `[1.2356 0 0.6322; 0 1.1266 0; 0.6322 0 8.2230] × 10^7`; units are not annotated in the source. The basis is `sphten-liouv` without approximation, and trajectory-level SSR is disabled.
- `powder` calls `@endor_mims` in the `esr` context using `rep_2ang_12800pts_sph`. The FID uses 128 points, sweep 3×10^8, τ = 100 ns and zero-filling to 512 points; `axis_units` is set to MHz.
- Before Fourier transformation, the mean is removed and the FID is apodised with the source setting `{'exp',6}`. The source then computes `fftshift(fft(fid,zerofill))` and plots the real spectrum against nuclear frequency in MHz with `plot_1d`. The function displays a figure and does not save the FID or spectrum.

## Requirements and scope

Run with MATLAB and Spinach, including `powder`, spin-system/basis construction, `endor_mims` and `apodisation`; plotting uses `plot_1d` and `kxlabel`. This is the ideal-hard-pulse Mims nitroxide workflow, not the long soft-pulse Davies calculation. The source provides no DOI.

Source: [examples/esr_sol_pulsed/endor_mims_nox_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/endor_mims_nox_powder.m).
