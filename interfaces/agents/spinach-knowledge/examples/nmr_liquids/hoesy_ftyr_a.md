# examples/nmr_liquids/hoesy_ftyr_a.m

- Signature: `hoesy_ftyr_a()`

## Purpose

(1H) -> (19F) HOESY spectrum of fluorotyrosine, with the magnetisation transfer direction picked so as to minimise the time that 19F spends in the transverse plane. This is the only way to run this sequence in proteins because aromatic 19F T2 is short. Calculation time: minutes

## Physical / mathematical content

This is the 1H-to-19F fluorotyrosine HOESY variant. The source loads the 3-fluorotyrosine DFT spin system and uses Redfield relaxation with a 10e-9 s correlation time and temperature 298 K, reflecting the stated large-protein model.

## Numerical / algorithmic content

The IK-2 sphten-liouv basis uses scalar-coupling connectivity and proximity level 3. A 0.5 s mixing time is simulated with 128 points per dimension, zero-filled to 512, followed by square-cosine apodisation and 2D Fourier transforms.

## Implementation structure

The script sets the field to 14.1 T and offsets [3000 -70000] Hz with sweeps [4000 2500] Hz. It observes spins `{'1H','19F'}`, decouples 19F in F1, and plots the real spectrum with negative display polarity.
