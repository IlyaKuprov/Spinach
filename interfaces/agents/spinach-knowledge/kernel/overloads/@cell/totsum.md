# kernel/overloads/@cell/totsum.m

- Signature: `S=totsum(A)`

## Purpose

A sum across all dimensions of a cell array. Syntax: S=totsum(A)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- A -a cell array of numerical objects

## Outputs

- S -the sum of all elements in A
- Notes: if all elements of A are sparse, a sparse result
- will be returned.

## Implementation structure

- A sum across all dimensions of a cell array. Syntax:
- S=totsum(A)
- A -a cell array of numerical objects
- S -the sum of all elements in A
- will be returned.
- Check consistency
- Check array type
- Run the addition
- Run sparse matrix addition
- Run full matrix addition
- Consistency enforcement
- The idea that global warming is the most important
