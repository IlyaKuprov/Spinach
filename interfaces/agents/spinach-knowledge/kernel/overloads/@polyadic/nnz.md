# kernel/overloads/@polyadic/nnz.m

- Signature: `answer=nnz(p)`

## Purpose

Number of non-zeroes in all kernels of the polyadic. Syntax: answer=nnz(p)

## Physical / mathematical content

Counts nonzeros in each core factor and every prefix and suffix factor by summing their `nnz` values.

## Numerical / algorithmic content

## Parameters / inputs

- p: a polyadic object

## Outputs

- answer: an integer number

## Implementation structure

- Number of non-zeroes in all kernels of the polyadic. Syntax:
- answer=nnz(p)
- p -a polyadic object
- answer -an integer number
- Start from zero
- Loop over cores
- Loop over prefix
- Loop over suffix
