# examples/relaxation_theory/maz_noesy_1.m

- Signature: `maz_noesy_1()`

## Purpose

Simulates a NOESY experiment for 15N-labelled methylaziridine, including scalar relaxation of the first kind from nitrogen-centre inversion modulating J-couplings. The source describes the effect in [the cited study](http://dx.doi.org/10.1002/ange.201410271) and gives a calculation time of minutes.

## Physical / mathematical content

The eight-spin model contains seven protons and one 15N nucleus at 11.75 T. Vacuum-DFT shielding tensors and scalar couplings are supplied; isotropic shifts are adjusted to the listed experimental values. Cartesian coordinates are in angstroms. The calculation combines Redfield relaxation and SRFK, with zero equilibrium, kite retention, a 25 ps correlation time, SRFK correlation parameters `[1.0 1e-3]`, and modulation depths 15 for couplings `(1,5)`, `(2,5)`, and `(3,5)`.

## Numerical / algorithmic content

The simulation uses an `sphten-liouv` basis with IK-2 approximation, scalar-coupling connectivity, and proximity level 3. Inter-spin and proximity cutoffs are 2.0 and 4.0, respectively, and Krylov propagation is disabled. `liquid` runs the `noesy` sequence with 2.0 s mixing time, 500 Hz offset, 1400 Hz sweeps in both dimensions, 256 points per dimension, and zero-filling to 1024 per dimension. Cosine apodisation is applied to both components of the FID before the two-dimensional Fourier transforms and plotting.

## Implementation structure

The script defines the eight isotopes, shielding matrices, experimental shifts, couplings, and coordinates, then configures the basis and both relaxation mechanisms. It sets the 1H pulse-sequence parameters and initial `Lz` state, runs `liquid(spin_system,@noesy,...)`, apodises the cosine and sine FIDs, transforms each component along F2, combines them as the states signal, transforms along F1, and plots the negative real spectrum.
