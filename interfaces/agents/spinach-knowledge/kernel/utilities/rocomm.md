# kernel/utilities/rocomm.m

- Signature: `C=rocomm(A)`

## Purpose

Right-ordered nested commutator [[[[A{1},A{2}],A{3}],A{4}],...] built from the user-supplied matrices. Syntax: C=rocomm(A)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A -a cell array of square matrices

## Outputs

- C -right-ordered nested commutator

## Implementation structure

- Right-ordered nested commutator [[[[A{1},A{2}],A{3}],A{4}],...]
- built from the user-supplied matrices. Syntax:
- C=rocomm(A)
- A -a cell array of square matrices
- C -right-ordered nested commutator
- Check consistency
- Nest the commutators
- Consistency enforcement
- It's not science I don't trust -it's the scientists.
- James Delingpole,
- a climate sceptic
