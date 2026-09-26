# kernel/utilities/kronm_new.m

- Signature: `M=kronm_new(Q,M)`

## Purpose

Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening Kronecker products. Syntax: M=kronm(Q,M)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- Q -cell array of Kronecker terms
- M -a vector or a matrix of appropriate dimension
- Output:
- M -a vector or a matrix of appropriate dimension

## Implementation structure

- Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening
- Kronecker products. Syntax:
- M=kronm(Q,M)
- Q - cell array of Kronecker terms
- M - a vector or a matrix of appropriate dimension
- Output:
- Check consistency
- Dimension statistics
- Row and column counts in Q
- Fold up implicit dimensions of M
- Run the products
- Contract each implicit dimension
