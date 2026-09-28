# kernel/utilities/distrib_dim.m

- Signature: `A=distrib_dim(A,dim)`

## Purpose

Distributes an array in the user-specified dimension for parallel processing using spmd. Syntax: A=distrib_dim(A,dim)

## Physical / mathematical content

- This function partitions a numerical array along a requested dimension; no physical model is specified.

## Numerical / algorithmic content

- Uses `spmd` and `codistributor1d(dim)` to partition `A` across workers along dimension `dim`, returning a distributed array.

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
