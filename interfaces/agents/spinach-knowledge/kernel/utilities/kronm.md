# kernel/utilities/kronm.m

- Signature: `x=kronm(Q,x)`

## Purpose

Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening Kronecker products. Syntax: x=kronm(Q,x)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- Q -cell array of Kronecker terms
- x -a vector or a matrix of appropriate dimension
- Output:
- x -a vector or a matrix of appropriate dimension

## Implementation structure

- Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening
- Kronecker products. Syntax:
- x=kronm(Q,x)
- Q - cell array of Kronecker terms
- x - a vector or a matrix of appropriate dimension
- Output:
- Check consistency
- Number of matrices in Q
- Number of columns in x
- Row and column counts in Q
- Dimension map for x
- Reshape into the map
