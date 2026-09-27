# examples/esr_sol_pulsed/endor_davies_nox_powder.m

- Signature: `endor_davies_nox_powder()`

## Purpose

Powder Davies ENDOR simulation for an electron–`14N` nitroxide system. The source models soft pulses with the Fokker–Planck formalism.

## Physical / mathematical content

- The calculation uses a spherical orientation grid and a `sphten-liouv` basis, with the electron Zeeman and electron–nitrogen coupling matrices defined in the source.

## Numerical / algorithmic content

- The example calls `powder(...,@endor_davies,...,'esr')` using the `rep_2ang_12800pts_sph` grid.
- It sweeps the nuclear pulse frequency from −80 to +80 MHz over 200 points and plots the real response.
