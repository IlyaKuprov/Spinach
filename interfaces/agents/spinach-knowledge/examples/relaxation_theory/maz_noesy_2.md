# examples/relaxation_theory/maz_noesy_2.m

- Signature: `maz_noesy_2()`

## Purpose

Simulates a methylaziridine NOESY spectrum with both scalar relaxation of the first kind, from nitrogen-centre inversion modulating J-couplings, and of the second kind, from rapid 14N quadrupolar relaxation. The source points to [the study describing this effect](http://dx.doi.org/10.1002/ange.201410271) and estimates minutes of calculation time.

## Physical / mathematical content

The eight-spin system contains seven protons and one 14N nucleus at 11.75 T. It uses vacuum-DFT shielding tensors, quadrupole and scalar couplings, and Cartesian coordinates in angstroms; isotropic shifts are assigned from experiment. Redfield, SRFK, and SRSK relaxation are enabled with zero equilibrium, kite retention, a 25 ps correlation time, and spin 4 as the SRSK source. SRFK uses correlation parameters `[1.0 1e-3]` and modulation depths of 15 for couplings `(1,5)`, `(2,5)`, and `(3,5)`.

## Numerical / algorithmic content

The calculation uses `sphten-liouv` with IK-2 approximation, scalar-coupling connectivity, and proximity level 4; inter-spin and proximity cutoffs are 2.0 and 4.0, and Krylov propagation is disabled. The NOESY sequence uses 2.0 s mixing, 500 Hz offset, 1400 Hz sweeps in both dimensions, 256 points per dimension, and 1024-point zero-filling in each dimension. Cosine apodisation and two-dimensional Fourier transforms produce the plotted spectrum.

## Implementation structure

The script specifies the eight nuclei, shielding and coupling tensors, shifts, and coordinates, configures the three relaxation contributions, and builds the spin system and basis. It runs `liquid` with `@noesy`, initial proton `Lz` state, and the NMR context. The cosine and sine FIDs are apodised and transformed along F2, combined into the states signal, transformed along F1, and plotted as the negative real spectrum.
