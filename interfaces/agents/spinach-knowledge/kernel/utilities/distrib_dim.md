# kernel/utilities/distrib_dim.m

- Signature: `A=distrib_dim(A,dim)`

## Purpose

Distributes an array in the user-specified dimension for parallel processing using spmd. Syntax: A=distrib_dim(A,dim)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- A -a numerical array
- dim -the distribution dimension
- Output:
- A -a distributed numerical array
- Mathworks, Inc.

## Implementation structure

- Distributes an array in the user-specified dimension
- for parallel processing using spmd. Syntax:
- A=distrib_dim(A,dim)
- A -a numerical array
- dim -the distribution dimension
- Output:
- A -a distributed numerical array
- Mathworks, Inc.
- Check consistency
- Get the size
- Set the stage
- Codistributor with default partitioning
