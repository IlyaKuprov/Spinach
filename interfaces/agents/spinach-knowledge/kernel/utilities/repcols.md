# kernel/utilities/repcols.m

## Purpose

Replicates specified columns of a matrix or cell array a specified number of times, returning a new array of the same type as the input.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/repcols.m>

## Behaviour

The function is called as `B=repcols(A,col_nums,rep_counts)`. It first runs a consistency check (`grumble`) on the inputs, then builds a replication map: a row vector of ones with one entry per column of `A`, where the entries at positions `col_nums` are set to the corresponding values of `rep_counts`. A column index vector is generated with `repelem(1:n,rep_map)`, where `n` is the number of columns of `A`, and the output is formed by indexing `A` with that column vector (`B=A(:,col_idx)`). Columns not listed in `col_nums` are kept exactly once; listed columns appear the number of times given by their replication count, in the order induced by the index vector.

The consistency check enforces:

- `A` must be numeric or a cell array, otherwise an error `'A must be numeric or a cell array.'` is thrown.
- `col_nums` must be a numeric vector of positive integers (no values below 1, no non-integers), otherwise an error `'col_nums must be a vector of positive integers.'` is thrown.
- `rep_counts` must be a numeric vector of positive integers with the same number of elements as `col_nums`, otherwise an error `'rep_counts must match col_nums in size and contain positive integers.'` is thrown.
- `col_nums` must not contain duplicates, otherwise an error `'col_nums must not contain duplicates.'` is thrown.
- `col_nums` indices must not exceed the number of columns of `A`, otherwise an error `'col_nums indices exceed matrix dimension.'` is thrown.

## Inputs and outputs

Inputs:

- `A` — a numeric matrix or a cell array.
- `col_nums` — vector of column indices to replicate.
- `rep_counts` — vector of positive integers specifying how many copies of each column to make; must match `col_nums` in size.

Output:

- `B` — same type as `A`, containing the replicated columns.

## References

- Spinach Wiki page for `repcols.m`: <https://spindynamics.org/wiki/index.php?title=repcols.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/repcols.m>
