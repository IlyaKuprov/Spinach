# kernel/utilities/keep_rank.m

## Purpose

Truncates the singular value decomposition of a matrix at a specified rank and reassembles the matrix, returning a low-rank approximation ([source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/keep_rank.m)).

## Behaviour

- Syntax: `A=keep_rank(A,nsvk)`.
- Runs a consistency check (`grumble`) on the inputs before processing.
- Converts the input to full storage with `full(A)` and computes the singular value decomposition `[U,S,V]=svd(full(A))`.
- Truncates the decomposition to the specified rank and rebuilds the matrix as `A=U(:,1:nsvk)*S(1:nsvk,1:nsvk)*V(:,1:nsvk)'`.
- The consistency check errors with `'A must be a matrix.'` if `A` is not numeric or if either dimension has size less than or equal to 1.
- The consistency check errors with `'nsvk must be a positive integer smaller than dim(A)'` if `nsvk` is not numeric, not real, not scalar, less than 1, non-integer (`mod(nsvk,1)~=0`), or greater than any dimension of `A` (`any(nsvk>size(A))`).

## Inputs and outputs

**Inputs**

- `A` — real or complex matrix; sparse inputs are converted to full.
- `nsvk` — number of singular values to keep.

**Outputs**

- `A` — filtered matrix, returned as full.

## References

- Source code: [kernel/utilities/keep_rank.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/keep_rank.m)
- Spinach Wiki: [keep_rank.m](https://spindynamics.org/wiki/index.php?title=keep_rank.m)
