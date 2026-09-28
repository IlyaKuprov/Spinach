# kernel/cache/st_product_table.m

- Signature: `[pt_left,pt_right]=st_product_table(nlevels)`

## Purpose

Structure coefficient tables for single transition operators.

## Physical / mathematical content

- The tables describe the left and right multiplicative action by `S{n}` on `S{m}`, as given in Eq 7.18 of the first edition of IK's book (normalisation is missing in the book, that's a typo).
- Numbering translation between single and double index is given by `kq2lin` and `lin2kq` functions.

## Numerical / algorithmic content

- The function obtains single transition operators using `sin_tran(nlevels)` and computes coefficients with `hdot(B{k},B{n}*B{m})` and `hdot(B{k},B{m}*B{n})`.
- These are expensive tables; a disk cache is used automatically. If a cache file exists, the function loads the tables instead of recomputing them. A failed attempt to save the cache produces a warning.

## Parameters / inputs

- `nlevels` - the number of energy levels in the system; it must be a positive integer scalar.

## Outputs

- `pt_left`, `pt_right` - structure coefficients in the following conventions:
  - `S{n}*S{m}=...+pt_left(n,m,k)*S{k}+...`
  - `S{m}*S{n}=...+pt_right(n,m,k)*S{k}+...`