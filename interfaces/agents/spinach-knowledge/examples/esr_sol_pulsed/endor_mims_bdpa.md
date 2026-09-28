# examples/esr_sol_pulsed/endor_mims_bdpa.m

- Signature: `endor_mims_bdpa()`

## Purpose

Mims ENDOR simulation on BDPA with ideal electron pulses, set to reproduce Figure 10 of [the cited paper](https://doi.org/10.1007/s00723-020-01269-z). The source estimates hours of calculation time and notes that a GPU can accelerate it.

## Physical / mathematical content

- The spin system contains one electron and two protons. The example uses a powder calculation for the Mims ENDOR sequence.

## Numerical / algorithmic content

- It sweeps 100 nuclear frequencies from 138 to 148 MHz with a 50 μs nuclear pulse, then plots absolute intensity against nuclear frequency in MHz.
