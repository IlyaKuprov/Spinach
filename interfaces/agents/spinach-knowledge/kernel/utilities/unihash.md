# kernel/utilities/unihash.m

- Signature: `A=unihash(A)`

## Purpose

Hash table based stable duplicate row eliminator, for use with large sparse matrices where Matlab's unique(...,'rows') is too slow. Syntax: A=unihash(A)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- A -a large and sparse matrix

## Outputs

- A -same matrix with duplicate
- rows deleted, keeping the
- first occurrence of each

## Implementation structure

- Hash table based stable duplicate row eliminator,
- for use with large sparse matrices where Matlab's
- unique(...,'rows') is too slow. Syntax:
- A=unihash(A)
- A -a large and sparse matrix
- A -same matrix with duplicate
- rows deleted, keeping the
- first occurrence of each
- Check consistency
- Build an MD5 hash table
- Redundant row index using a hash table
- Elimination
