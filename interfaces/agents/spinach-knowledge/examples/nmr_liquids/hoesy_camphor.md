# examples/nmr_liquids/hoesy_camphor.m

- Signature: `hoesy_camphor()`

## Purpose

13C{1H} HOESY spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed with DFT. Calculation time: minutes

## Physical / mathematical content

The example simulates a 13C-detected, 1H-decoupled HOESY experiment. It obtains the camphor spin system and interactions from the supplied DFT output, applies Redfield relaxation, and sums the simulated signal over 13C isotopomers.

## Numerical / algorithmic content

It uses an IK-2 sphten-liouv basis with scalar-coupling connectivity and proximity level 3. Each isotopomer is simulated with a 0.5 s mixing time; cosine apodisation and two Fourier transforms produce the 2D spectrum.

## Implementation structure

The script reads camphor data via `gparse`/`g2spinach`, sets a 14.1 T field and Redfield parameters (correlation time 50e-12 s, 298 K), then dilutes the system for 13C and runs a `parfor` loop. Acquisition uses sweeps [1800 9000] Hz, offsets [900 4500] Hz, 128 points per dimension and 512-point zero filling. The result is plotted with positive display polarity.
