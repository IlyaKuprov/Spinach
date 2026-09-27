# examples/esr_sol_pulsed/eseem_nitroxide_crystal.m

- Signature: `eseem_nitroxide_crystal()`

## Purpose

Two-pulse X-band ESEEM spectrum of a nitroxide radical at a specific orientation relative to the lab frame. Ideal pulses are assumed. Calculation time: seconds.

## Physical / mathematical content

- The spin system contains an electron and 14N at 0.33 T, with electron Zeeman and electron–nitrogen coupling matrices specified directly in the source.
- The crystal orientation is `[pi/5 pi/4 pi/3]`; no powder averaging is performed.

## Numerical / algorithmic content

- The simulation uses 1024 points at a 12.5 ns timestep. After mean subtraction and Kaiser apodisation with parameter 6, it applies an FFT with 4096-point zero filling and `fftshift`.
- The frequency axis uses an interpulse-delay increment of half the timestep.

## Implementation structure

- Create the spin system in the `sphten-liouv` basis without approximation, then call `crystal` with `@eseem` in the `esr` context.
- Plot the real time-domain signal and the magnitude spectrum.
