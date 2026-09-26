# kernel/utilities/sec2kite.m

- Signature: `R=sec2kite(spin_system,R)`

## Purpose

Converts a secular relaxation superoperator into the Redfield kite form by dropping all non-longitudinal cross-relaxation pro- cesses. Useful when the relaxation superoperator is huge, but TROSY-like effects are negligible. Syntax: R=sec2kite(spin_system,R)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

## Parameters / inputs

- R -relaxation superoperator

## Outputs

- R -relaxation superoperator

## Implementation structure

- Converts a secular relaxation superoperator into the Redfield
- kite form by dropping all non-longitudinal cross-relaxation pro-
- cesses. Useful when the relaxation superoperator is huge, but
- TROSY-like effects are negligible. Syntax:
- R=sec2kite(spin_system,R)
- R -relaxation superoperator
- Check consistency
- Get nonzero count
- Compile the index of all longitudinal product states in the basis
- Convert R to XYZ format
- Zero all rates except self-relaxation and longitudinal cross-relaxation terms
- Recompose the relaxation superoperator and get nonzero count
