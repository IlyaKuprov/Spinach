# kernel/utilities/unihash.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/unihash.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/unihash.m)

## Purpose

`unihash` is a hash-table-based stable duplicate row eliminator, intended for large sparse matrices where MATLAB's `unique(...,'rows')` is too slow. It removes duplicate rows from a matrix while keeping the first occurrence of each row.

## Behaviour

- Syntax: `A=unihash(A)`.
- The function first validates its input via an internal consistency check (`grumble`), which errors with `'A must be a numeric matrix.'` if the input is not numeric or not a matrix.
- An MD5 hash table is built as a character array of blanks, `repmat(' ',[size(A,1) 32])`, i.e. one 32-character row per row of `A`.
- A `parfor` loop over `k=1:size(A,1)` fills each row of the hash table with `md5_hash(A(k,:))`, so row hashing is parallelised.
- Redundant row indices are found with `[~,idx]=unique(hash_table,'rows','stable')`, which preserves the order of first occurrences.
- The elimination step returns `A=A(idx,:)`, deleting duplicate rows while keeping the first occurrence of each.

## Inputs and outputs

**Input:**

- `A` — a large and sparse matrix; must be numeric and a matrix.

**Output:**

- `A` — the same matrix with duplicate rows deleted, keeping the first occurrence of each.

## References

- Spinach Dynamics Wiki: [unihash.m](https://spindynamics.org/wiki/index.php?title=unihash.m)
