# kernel/utilities/swizzle.m

## Purpose

Flattens nested index lists into an array of tuples in random order, which is useful for flattening nested loops for parallel processing.

## Behaviour

- Syntax: `tuples=swizzle(index_arrays)`.
- The function first validates its input via an internal consistency check (`grumble`):
  - Errors with `index_arrays must be a cell array of row vectors.` if the input is not a cell array.
  - Errors with `elements of index_arrays must be row vectors of positive integers.` if any element is not real, not a row vector, contains non-integer values, or contains values less than 1.
- The tuples are built by Kronecker-style expansion: the first index vector initialises the column, and each subsequent vector appends a new column formed by `kron` of ones and the vector against the accumulated tuples.
- After construction, the rows of the tuple matrix are randomly permuted with `randperm`, so the tuples are returned in random order with tuples listed as rows.

## Inputs and outputs

**Inputs**

- `index_arrays` — a cell array of row vectors (each element must be a row vector of positive integers).

**Outputs**

- `tuples` — a matrix of tuples in random order, with tuples listed as rows.

## References

- Source: [kernel/utilities/swizzle.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/swizzle.m)
- Wiki: <https://spindynamics.org/wiki/index.php?title=swizzle.m>
