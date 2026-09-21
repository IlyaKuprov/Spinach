# kernel/utilities/rlx_split.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_split.m`
- Signature: `[R1,R2,Rm]=rlx_split(spin_system,R)`
- Total lines: 67

## Purpose

Splits a relaxation superoperator into longitudinal, trans- verse and mixed components. Syntax: [R1,R2,Rm]=rlx_split(spin_system,R)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `lin2lm()`, `logical()`, `any()`, `isfield()`, `strcmp()`, `ismatrix()`.
