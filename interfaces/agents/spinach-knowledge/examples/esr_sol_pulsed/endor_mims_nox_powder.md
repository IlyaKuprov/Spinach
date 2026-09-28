# examples/esr_sol_pulsed/endor_mims_nox_powder.m

- Signature: `endor_mims_nox_powder()`

## Purpose

Mims ENDOR simulation for a nitroxide radical powder. Ideal hard pulses are assumed. Calculation time: seconds.

## Physical / mathematical content

- The spin system contains an electron and 14N at 3.5 T. The electron g matrix has diagonal values 2.01045, 2.00641, and 2.00211; the source specifies an electron–nitrogen coupling matrix.
- Powder averaging uses `rep_2ang_12800pts_sph`.

## Numerical / algorithmic content

- The sequence uses 128 points, a 3e8 sweep parameter, a 100 ns interpulse delay, and 512-point zero filling.
- The signal is mean-subtracted, exponentially apodised with parameter 6, then Fourier transformed and FFT-shifted.

## Implementation structure

- Create the spin system in the `sphten-liouv` basis without approximation, with trajectory-level SSR disabled.
- Run `powder` with `@endor_mims` in the `esr` context.
- Plot the real spectrum against nuclear frequency in MHz.
