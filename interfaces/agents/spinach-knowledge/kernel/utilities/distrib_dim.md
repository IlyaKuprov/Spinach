# kernel/utilities/distrib_dim.m

## Purpose

Distributes a numerical array along a user-specified dimension for parallel processing using `spmd`, returning a distributed array suitable for parallel computation in Spinach.

## Behaviour

- The function first validates its inputs via an internal consistency check (`grumble`):
  - `A` must be numeric.
  - `dim` must be a real, scalar, positive integer.
  - `dim` must not exceed the number of dimensions of `A`.
- The size of `A` is obtained with `size`.
- Inside an `spmd` block, a 1-D codistributor is created with default (unset) partitioning: `codistributor1d(dim, codistributor1d.unsetPartition, size_A)`. A `LocalParts` variable is initialised to `1` on each worker.
- After the `spmd` block, the codistributor is taken from the first lab (`CoD{1}`), and partition limits are computed as `[0 cumsum(CoD.Partition)]`.
- A dimension-agnostic cell array of index colons (`repmat({':'},1,ndims(A))`) is built; the `dim`-th entry is replaced, for each local part `n`, by the index range `partLimits(n)+1:partLimits(n+1)`, and the corresponding slice of `A` is extracted into `LocalParts{n}`.
- The output is assembled with `A = distributed(LocalParts, dim)`, producing a distributed array whose local parts correspond to the computed partition along `dim`.

## Inputs and outputs

**Syntax:** `A = distrib_dim(A, dim)`

- **A (input)** — a numerical array to be distributed.
- **dim (input)** — the distribution dimension; a positive integer not exceeding `ndims(A)`.
- **A (output)** — a distributed numerical array partitioned along `dim`.

## References

- Source: [kernel/utilities/distrib_dim.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/distrib_dim.m)
- Spinach Wiki: [distrib_dim.m](https://spindynamics.org/wiki/index.php?title=distrib_dim.m)
- Attributed in the source header to Mathworks, Inc. and ilya.kuprov@weizmann.ac.il.
