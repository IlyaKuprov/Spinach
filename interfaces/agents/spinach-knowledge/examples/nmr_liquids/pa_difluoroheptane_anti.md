# examples/nmr_liquids/pa_difluoroheptane_anti.m

- Signature: `pa_difluoroheptane_anti()`

## Purpose

Pulse-acquire 1H NMR spectrum of anti-3,5-difluoroheptane. The example manually specifies a basis by merging Lie algebras for selected spin fragments, followed by symmetry factorisation and conservation-law screening. See [the cited paper](https://doi.org/doi/10.1021/acs.joc.4c00670) for further information. Calculation time: minutes, faster with a GPU.

## Physical / mathematical content

- The manually specified 23-spin model contains carbon-12, proton and fluorine-19 spins, with chemical shifts and scalar couplings assigned explicitly. The basis partitions the spins into three fragments and applies S3 symmetry to the two proton triplets; the longitudinal projection is specified for fluorine-19.

## Numerical / algorithmic content

- The basis uses spherical-tensor Liouville formalism, IK-0 approximation and interaction level 1. Manual fragment membership, symmetry factorisation and projection selection define the reduced basis.
- GPU acceleration is noted as useful but remains commented out in the source; ZTE is explicitly disabled. Acquisition uses 4096 points, zero-filling to 16536, a 2500 Hz sweep and 1400 Hz offset, followed by 5 Hz exponential apodisation.

## Implementation structure

- Define the field (11.7464 T), 23 isotopes, chemical shifts and scalar couplings, then construct the manually partitioned, symmetry-factorised basis.
- Set proton pulse-acquire initial state and coil, simulate, apodise, Fourier transform and plot the real spectrum with the frequency axis inverted.
