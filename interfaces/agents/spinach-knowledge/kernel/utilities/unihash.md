# kernel/utilities/unihash.m

- Signature: `A=unihash(A)`

## Purpose

Stable duplicate-row eliminator based on a hash table, intended for large sparse matrices where MATLAB's `unique(...,'rows')` is too slow.

## Physical / mathematical content

This utility removes duplicate rows and retains the first occurrence of each row.

## Numerical / algorithmic content

The function computes an MD5 hash for each row in a `parfor` loop, uses stable row-wise uniqueness of the hashes to identify redundant rows, and returns the retained rows in their original order.

## Parameters / inputs

- `A` — large sparse matrix.

## Outputs

- `A` — the same matrix with duplicate rows deleted, keeping the first occurrence of each.

## Implementation structure

- Checks the input.
- Builds an MD5 hash table in parallel.
- Finds redundant rows using stable hash-table uniqueness and removes them.
- Source documentation: <https://spindynamics.org/wiki/index.php?title=unihash.m>
