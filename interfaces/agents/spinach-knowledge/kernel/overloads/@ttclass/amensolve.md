# kernel/overloads/@ttclass/amensolve.m

- Signature: `x=amensolve(A,y,tol,opts,x0)`

## Purpose

Approximately solve the tensor-train linear system `A*x=y` with the alternating minimal energy (AMEn) iteration.

## Physical / mathematical content

The matrix, right-hand side, initial guess, and solution are represented as tensor trains. The routine requires a square matrix and a vector right-hand side with matching mode sizes; matrix and right-hand side must each contain one train and be shrunk before solving.

## Parameters / inputs

- `A` - `ttclass` representing a square matrix.
- `y` - `ttclass` right-hand-side vector with matching mode sizes.
- `tol` - finite, nonnegative real scalar relative tolerance; the source suggests `1e-6` as a starting value.
- `opts` - optional options structure (an empty value selects defaults):
  - `nswp` - maximum AMEn sweeps (default `20`).
  - `init_guess_rank` - rank used for the random initial guess (default `2`).
  - `enrichment_rank` - residual/enrichment rank (default `4`; set to zero to disable enrichment).
  - `resid_damp` - local accuracy damping factor (default `2`).
  - `rmax` - maximum solution TT rank (default `Inf`).
  - `max_full_size` - local problem size below which a direct solve is used (default `500`).
  - `local_iters` - maximum BiCGSTAB iterations for local problems (default `100`).
  - `sv_floor` - absolute singular-value floor used during solution-block truncation (default `1e-10`). Because it is absolute, tune it to the scale of the orthogonalised solution blocks, which follows `norm(y)/norm(A)`; lower it when the solution is small. Any extra truncation from raising the floor is capped by and charged to the block's truncation budget; use `tol` to request a coarser returned solution beyond that cap.
  - `verb` - verbosity level: silent (0), sweep (1), or full (2); default `1`.
- `x0` - optional `ttclass` initial-guess vector with matching mode sizes; if omitted, a random initial guess is generated.

## Outputs

- `x` - `ttclass` approximate solution. The source documents the target `|x-X| < tol*|X|` in Frobenius norm for exact solution `X`; the returned train records a tolerance.

## Implementation structure

The AMEn sweeps solve local problems directly below `max_full_size`, otherwise with BiCGSTAB. The solution blocks are truncated by SVD subject to the tolerance, singular-value floor, and rank limit. When enabled, the projected residual enriches the train; sweeps stop when the iterate change is below `tol` or `nswp` is reached.
