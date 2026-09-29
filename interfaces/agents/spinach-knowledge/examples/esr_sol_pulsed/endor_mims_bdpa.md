# examples/esr_sol_pulsed/endor_mims_bdpa.m

Call `endor_mims_bdpa()` with no arguments. This powder Mims ENDOR example models BDPA with ideal electron pulses and is intended to reproduce Figure 10 of [the cited paper](https://doi.org/10.1007/s00723-020-01269-z). The source estimates hours of calculation and notes it is much faster on a GPU.

## Spin system and sequence

- The system contains one electron and two 1H nuclei at 3.35 T. The diagonal electron g values are 2.00263, 2.00260 and 2.00257. The two electron–proton coupling matrices are diagonal: `[7.70, 5.30, 2.00] × 10^6` and `[1.00, 1.00, 1.26] × 10^6`; matrix units are not labelled in the source.
- The source uses a `sphten-liouv` basis without approximation, `t1_t2` relaxation, zero equilibrium and diagonal relaxation retention. Its `r1_rates` are `{1e3, 1e4, 1e4}` and `r2_rates` are `{1e1, 1e2, 1e2}` (rate units are not stated in the file).
- `powder` runs `@endor_mims_ideal` in the `esr` context, with `rep_2ang_400pts_sph`. The sequence parameters select the electron, specify nuclei 2 and 3, and set τ = 200 ns. The nuclear π-pulse duration is 50 μs; its RF field is set by `−π/(n_dur × spin('1H'))`. The nuclear-frequency sweep has 100 points from 138 to 148 MHz, with nuclear pulse rank 2.
- The displayed trace is `abs(answer)` against laboratory-frame nuclear frequency in MHz; its y-axis is absolute intensity in arbitrary units. The example produces a figure, not a saved data file.

## Requirements and scope

Run under MATLAB with Spinach, including its `powder`, spin-system/basis construction and `endor_mims_ideal` sequence code; the plotting code uses `kfigure` and `kgrid`. The ideal-electron-pulse assumption is part of this example, not a finite electron-pulse simulation.

Source: [examples/esr_sol_pulsed/endor_mims_bdpa.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/endor_mims_bdpa.m).
