# examples/nmr_liquids/hsqc_sucrose.m

- Signature: `hsqc_sucrose()`

## Purpose

HSQC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time: seconds

## Physical / mathematical content

Magnetic parameters are generated from the sucrose vacuum-DFT output. Before simulation, the script replaces selected isotropic shielding values with experimental shifts; the HSQC is then accumulated over 13C isotopomers.

## Numerical / algorithmic content

The IK-2 sphten-liouv basis uses scalar-coupling connectivity and proximity level 1. The sequence uses J=140 Hz, square apodisation, 128 points per dimension and 512-point zero filling; simulations are run in parallel and combined using the States signal before the 2D Fourier transform.

## Implementation structure

The field is 5.9 T; sweeps are [3350 950] Hz and offsets [5000 1100] Hz. The script updates shielding entries for spin numbers [1:19 24:30] using the listed `new_shifts` array, dilutes for 13C, and plots the real spectrum with positive polarity.
