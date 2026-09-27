# examples/esr_sol_pulsed/endor_davies_nox_crystal.m

- Signature: `endor_davies_nox_crystal()`

## Purpose

Two-stage single-orientation nitroxide Davies ENDOR example. Soft pulses are simulated with the Fokker–Planck formalism.

## Physical / mathematical content

- The first stage computes a crystal-orientation pulse-acquire ESR spectrum; the second runs the Davies ENDOR sequence.
- For the ENDOR stage, the electron-frequency offsets are +93, +10, and −73 MHz, and the nuclear-frequency sweep spans −200 to +200 MHz.

## Numerical / algorithmic content

- The ESR FID is apodised, Fourier transformed with 512-point zero filling, and plotted. The ENDOR simulation uses the source-defined frequency sweep and crystal-orientation setup.
