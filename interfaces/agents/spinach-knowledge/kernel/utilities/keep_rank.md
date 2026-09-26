# kernel/utilities/keep_rank.m

- Signature: `A=keep_rank(A,nsvk)`

## Purpose

Truncates the singular value decomposition at the specified rank and reassembles the matrix. Syntax: A=keep_rank(A,rank)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A -real or complex matrix, will be
- converted to full if a sparse
- matrix is received
- nsvk -number of singular values to keep

## Outputs

- A -filtered matrix, returned as full

## Implementation structure

- Truncates the singular value decomposition at the specified rank
- and reassembles the matrix. Syntax:
- A=keep_rank(A,rank)
- A - real or complex matrix, will be
- converted to full if a sparse
- matrix is received
- nsvk - number of singular values to keep
- A - filtered matrix, returned as full
- Check consistency
- Run singular value decomposition
- Truncate to the specified rank and rebuild
- Consistency enforcement
