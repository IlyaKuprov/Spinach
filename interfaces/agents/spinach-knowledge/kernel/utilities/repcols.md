# kernel/utilities/repcols.m

- Signature: `B=repcols(A,col_nums,rep_counts)`

## Purpose

Replicates specified columns of a matrix or cell array a specified number of times. Syntax: B=repcols(A,col_nums,rep_counts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A -a numeric matrix or a cell array
- col_nums -vector of column indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each column to make
- Output:
- B -same type as A

## Implementation structure

- Replicates specified columns of a matrix or cell array a
- specified number of times. Syntax:
- B=repcols(A,col_nums,rep_counts)
- A -a numeric matrix or a cell array
- col_nums -vector of column indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each column to make
- Output:
- B -same type as A
- Check consistency
- Replication counts for every column
- Build column index vector
