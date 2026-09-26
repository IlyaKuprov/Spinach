# kernel/utilities/rlx_split.m

- Signature: `[R1,R2,Rm]=rlx_split(spin_system,R)`

## Purpose

Splits a relaxation superoperator into longitudinal, trans- verse and mixed components. Syntax: [R1,R2,Rm]=rlx_split(spin_system,R)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- R -a relaxation superoperator in sphten-liouv
- formalism

## Outputs

- R1 -the part of R acting on purely longitudinal
- single-spin states
- R2 -the part of R acting on purely transverse
- single-spin states
- Rm -the rest of R

## Implementation structure

- Splits a relaxation superoperator into longitudinal, trans-
- verse and mixed components. Syntax:
- [R1,R2,Rm]=rlx_split(spin_system,R)
- R -a relaxation superoperator in sphten-liouv
- formalism
- R1 -the part of R acting on purely longitudinal
- single-spin states
- R2 -the part of R acting on purely transverse
- Rm -the rest of R
- Check consistency
- Interpret the basis
- Index single-spin orders
