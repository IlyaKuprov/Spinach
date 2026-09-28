# kernel/overloads/@ttclass/nnz.m

- Signature: `answer=nnz(ttrain)`

## Purpose

Counts nonzero entries across the cores of a tensor train.

## Physical / mathematical content

For each core, MATLAB `nnz` counts its nonzero entries; the function sums those counts over all cores. It does not expand the represented tensor.

## Numerical / algorithmic content

The implementation applies `nnz` to each cell in `ttrain.cores`, then sums the resulting counts.

## Parameters / inputs

- ttrain - tensor train object

## Outputs

- answer - number of nonzero elements across all tensor train cores

## Implementation structure

- Apply `nnz` to every core with `cellfun`.
- Sum the per-core counts.
