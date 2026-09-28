# examples/nmr_solids/dor_powder_nav_fplanck_time.m

- Signature: `dor_powder_nav_fplanck_time()`

## Purpose

Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calculation includes the second-order quadrupolar shift and the third-order lineshape. Time-domain detection. Note: slower spinning rates and larger NQIs require larger ranks and spherical grids. At the moment the spinning frequencies are set artificially too high to reduce the simulation time in this example. Calculation time: seconds

## Physical / mathematical content

This 14N double-angle-spinning example models the quadrupolar interaction with `eeqq2nqi(3.21e6,0.27,1,[0 0 0])` at 14.1 T. The source describes the target as the 14N nucleus of N-acetylvaline and notes second-order quadrupolar shift and third-order lineshape contributions.

## Numerical / algorithmic content

The source calls `doublerot` with `acquire` in the lab frame, using the 1D Fokker–Planck treatment, outer/inner rates of 1 and 5 MHz and ranks 7 and 4. It uses the `rep_2ang_100pts_oct` grid and disables `trajlevel`; diagonal damping is 2 kHz. The time-domain signal has 256 points over a 100 kHz sweep and is zero-filled to 1024 before Fourier transformation.

## Implementation structure

Builds the single-spin quadrupolar system and basis, sets DOR and acquisition parameters, acquires the time-domain signal, Fourier transforms it, and plots the real spectrum.
