# kernel/overloads/@cell/totsum.m

- Signature: `S=totsum(A)`

## Purpose

Elementwise sum of the numeric arrays in a cell array. Syntax: S=totsum(A)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- A -a cell array of numerical objects

## Outputs

- S -the elementwise sum of the numeric arrays in A
- Notes: if all elements of A are sparse, a sparse result
- will be returned.

## Implementation structure

- Check consistency
- Check array type
- Collect nonzeros and construct a sparse result when all entries are sparse
- Run full matrix addition otherwise
- Enforce numeric-array consistency
