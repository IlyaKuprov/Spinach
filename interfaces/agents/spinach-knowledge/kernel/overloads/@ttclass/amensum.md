# kernel/overloads/@ttclass/amensum.m

- Signature: `y=amensum(x,tol,opts)`

## Purpose

Compress a sum of buffered rank-one tensor trains into one tensor train using an AMEn iteration.

## Physical / mathematical content

The input represents a CP-format sum: each core has rank one across the buffered trains. The routine currently handles this case; buffered sums of general tensor trains are not implemented.

## Parameters / inputs

- `x` - `ttclass` containing buffered rank-one tensor trains.
- `tol` - relative tolerance used for convergence and SVD truncation (the source suggests `1e-10`).
- `opts` - optional options structure:
  - `max_swp` - maximum iteration count (default `100`).
  - `init_guess_rank` - rank of the random initial guess (default `2`).
  - `enrichment_rank` - rank of residual enrichment (default `4`; zero disables enrichment).
  - `verb` - verbosity switch (default `0`).

## Outputs

- `y` - one `ttclass` tensor train, with the documented target `|x-y| < tol*|x|` in Frobenius norm.

## Implementation structure

The routine initializes a random tensor train, updates its cores by alternating local projections of the buffered rank-one terms, and truncates intermediate cores by SVD. Optional residual enrichment augments the approximation. Iteration stops when the largest relative core change is below `tol` or `max_swp` is reached. Non-rank-one buffered input currently raises an error.
