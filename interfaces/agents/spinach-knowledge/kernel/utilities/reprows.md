# kernel/utilities/reprows.m

- Signature: `B=reprows(A,row_nums,rep_counts)`

## Purpose

Replicates specified rows of a matrix or cell array a specified number of times. Syntax: B=reprows(A,row_nums,rep_counts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A -a numeric matrix or a cell array
- row_nums -vector of row indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each row to make

## Outputs

- B -same type as A

## Implementation structure

- Replicates specified rows of a matrix or cell array a
- specified number of times. Syntax:
- B=reprows(A,row_nums,rep_counts)
- A -a numeric matrix or a cell array
- row_nums -vector of row indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each row to make
- B -same type as A
- Check consistency
- Replication counts for every row
- Build row index vector
- Extract and replicate
