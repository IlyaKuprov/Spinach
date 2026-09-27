# examples/nmr_liquids/hoesy_ftyr_b.m

- Signature: `hoesy_ftyr_b()`

## Purpose

(19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not the right way to run this sequence in proteins because aromatic 19F T2 is short, but 19F is being phase-encoded. Calculation time: minutes

## Physical / mathematical content

This is the reverse 19F-to-1H fluorotyrosine HOESY variant. The source loads the 3-fluorotyrosine DFT spin system and uses Redfield relaxation with a 10e-9 s correlation time and temperature 298 K, with the code comment identifying the correlation time as appropriate to a large protein.

## Numerical / algorithmic content

The IK-2 sphten-liouv basis uses scalar-coupling connectivity and proximity level 3. A 0.5 s mixing time is simulated with 128 points per dimension, zero-filled to 512, then square-cosine apodisation and 2D Fourier transforms are applied.

## Implementation structure

The script sets a 14.1 T field, sweeps [2500 4000] Hz and offsets [-70000 3000] Hz. It observes spins `{'19F','1H'}`, decouples 1H in F1, and plots the real spectrum with negative display polarity.
