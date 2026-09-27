# examples/nmr_solids/dor_powder_nav_fplanck_freq.m

- Signature: `dor_powder_nav_fplanck_freq()`

## Purpose

Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calculation includes the second-order quadrupolar shift and the third-order lineshape. Frequency-domain detection within the user-specified frequency interval. Note: slower spinning rates and larger NQIs require larger ranks and spherical grids. At the moment the spinning frequencies are set artificially too high to reduce the simulation time in this example. Calculation time: minutes

## Physical / mathematical content

This 14N double-angle-spinning example models the quadrupolar interaction with `eeqq2nqi(3.21e6,0.27,1,[0 0 0])` at 14.1 T. The source identifies the target as the 14N nucleus of N-acetylvaline and describes the spectrum as including second-order quadrupolar shift and third-order lineshape contributions.

## Numerical / algorithmic content

`doublerot` runs the `slowpass` sequence in the lab frame with the 1D Fokker–Planck treatment, outer/inner rates of 1 and 5 MHz, and ranks 7 and 4. The octahedral spherical grid is `rep_2ang_100pts_oct`; the calculated frequency-domain spectrum spans −50 to +50 kHz at 1024 points. The example disables `trajlevel` and includes diagonal damping at 2 kHz.

## Implementation structure

Creates the single-spin quadrupolar system and basis, sets the DOR rates, axes, ranks, grid and spectral interval, calculates the frequency-domain spectrum, and plots its real part.
