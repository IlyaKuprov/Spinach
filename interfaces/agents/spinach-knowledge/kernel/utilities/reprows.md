# kernel/utilities/reprows.m

## Purpose

Replicates specified rows of a matrix or cell array a specified number of times, as documented in the file header.

## Behaviour

The function first runs a consistency check (`grumble`) on the inputs. It then builds a replication map: a row vector of ones of length equal to the number of rows of `A`, in which the entries at positions given by `row_nums` are replaced by the corresponding values of `rep_counts`. A row index vector is generated with `repelem(1:n,rep_map)`, and the output is formed by indexing `A` with that vector across all columns (`B=A(row_idx,:)`).

The consistency check enforces:

- `A` must be numeric or a cell array, otherwise an error `'A must be numeric or a cell array.'` is raised.
- `row_nums` must be a real, finite numeric vector of positive integers, otherwise an error `'row_nums must be a vector of positive integers.'` is raised.
- `rep_counts` must be a real, finite numeric vector of positive integers with the same number of elements as `row_nums`, otherwise an error `'rep_counts must match row_nums in size and contain positive integers.'` is raised.
- `row_nums` must not contain duplicates, otherwise an error `'row_nums must not contain duplicates.'` is raised.
- `row_nums` entries must not exceed the number of rows of `A`, otherwise an error `'row_nums indices exceed matrix dimension.'` is raised.

## Inputs and outputs

Syntax: `B=reprows(A,row_nums,rep_counts)`

- `A` — a numeric matrix or a cell array.
- `row_nums` — vector of row indices to replicate.
- `rep_counts` — vector of positive integers specifying how many copies of each row to make.
- `B` — output, same type as `A`.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/reprows.m>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=reprows.m>
