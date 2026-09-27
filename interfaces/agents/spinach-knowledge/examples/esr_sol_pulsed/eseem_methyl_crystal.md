# examples/esr_sol_pulsed/eseem_methyl_crystal.m

- Signature: `eseem_methyl_crystal()`

## Purpose

Two-pulse X-band ESEEM spectrum of a methyl radical at a specific orientation relative to the lab frame. Magnetic parameters are imported from a vacuum-DFT calculation. Ideal pulses are assumed. Calculation time: seconds.

## Physical / mathematical content

- Spin-system properties are imported from `../standard_systems/methyl.log`, mapping the electron and hydrogen to `E` and `1H`.
- The magnetic field is 0.33 T and the crystal orientation is `[pi/5 pi/4 pi/3]`; no powder averaging is performed.

## Numerical / algorithmic content

- The simulation uses 512 points at a 10 ns timestep. The mean-subtracted signal receives Kaiser apodisation with parameter 6 before an FFT with 4096-point zero filling and `fftshift`.

## Implementation structure

- Create the spin system in the `sphten-liouv` basis without approximation, then call `crystal` with `@eseem` in the `esr` context.
- Plot the real time-domain signal and the magnitude spectrum.
