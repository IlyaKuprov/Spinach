# examples/spin_chemistry/singlet_yield_nz_2.m

- Signature: `singlet_yield_nz_2()`

## Purpose

Compare Redfield and lifetime-shifted Nakajima-Zwanzig (NZ) predictions for the field-dependent decay rate of a micelle-confined, triplet-born benzophenone ketyl/alkyl radical pair. The source describes high-field decay as relaxation-controlled: anisotropic hyperfine modulation transfers population from the T+/- states into the reactive S/T0 subspace. The reported decay rate is the slowest eigenmode of the pair Liouvillian.

## Setup

The calculation uses fields `0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.34 T`, contact recombination rate `1e9 Hz`, and rotational correlation time `0.7 ns`. The source comments place its scalar lifetime-shift parameter near `0.35` and identify the SDS supercage systems of Sakaguchi, Hayashi, and Nagakura as the experimental context. It calculates both theories at every field and plots the two decay-rate curves.

## Reference

- [DOI: 10.1246/bcsj.57.322](https://doi.org/10.1246/bcsj.57.322)
