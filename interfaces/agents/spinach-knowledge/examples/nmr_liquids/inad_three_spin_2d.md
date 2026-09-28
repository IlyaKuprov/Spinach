# examples/nmr_liquids/inad_three_spin_2d.m

- Signature: `inad_three_spin_2d()`

## Purpose

Example of a 2D-INADEQUATE spectrum of a generic three-spin system. Calculation time: seconds. Contributors named in the source: Theresa Hune and Christian Griesinger.

## Implementation

The three 13C spins are simulated at 16.44 T (700 MHz) with J(1,2) = 20 Hz, J(1,3) = 60 Hz, and no coupling between spins 2 and 3. The code generates pair-labelled 13C isotopomers and simulates a two-dimensional INADEQUATE experiment. It apodises both quadrature components, forms the States signal, Fourier transforms both dimensions, and plots the resulting spectrum.
