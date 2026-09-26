# kernel/utilities/killcross.m

- Signature: `M=killcross(M,f1idx,f2idx)`

## Purpose

Zeroes the specified rows and columns of a matrix. Syntax: M=killcross(M,f1idx,f2idx)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- M -a matrix
- f1idx -numbers of the columns that
- should be zeroed
- f2idx -numbers of the rows that
- should be zeroed

## Outputs

- M -a matrix

## Implementation structure

- Zeroes the specified rows and columns of a matrix. Syntax:
- M=killcross(M,f1idx,f2idx)
- M -a matrix
- f1idx -numbers of the columns that
- should be zeroed
- f2idx -numbers of the rows that
- Check consistency
- Wipe the indices
- Consistency enforcement
- A narcissist is someone better-looking than you are.
- Gore Vidal
