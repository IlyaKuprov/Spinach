# examples/extremes/fluoroisooctane.m

- Signature: `fluoroisooctane()`

## Purpose

A deliberately adversarial example from Art Bochevarov at Schrödinger, Inc. In this case, IK-2 approximation in Liouville space generates an exceedingly large basis set; the calculation must instead be performed in Hilbert space with permutation symmetry factorisation. Calculation time: hours.

## Physical / mathematical content

- The target observable is the 1H NMR spectrum of a highly coupled fluoroisooctane spin system.
- The calculated observable is the proton free-induction decay and its Fourier-transformed NMR spectrum.

## Numerical / algorithmic content

- The IK-2 Liouville basis is impractically large, so the script uses Hilbert-space formalism with three S3 permutation-symmetry blocks before calculating and Fourier-transforming the proton FID.

## Implementation structure

- A deliberately adversarial example from Art Bochevarov at
- Schodinger Inc. In this case, IK-2 approximation in Liou-
- ville space generates an exceedingly large basis set; the
- calculation must instead be performed in Hilbert space
- with permutation symmetry factorisation.
- Calculation time: hours.
- Magnet induction
- Isotopes
- Chemical shifts
- Larger J-couplings
- Smaller J-couplings, tert-butyl
- Smaller J-couplings, isopropyl
