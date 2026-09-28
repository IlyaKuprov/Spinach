# kernel/overloads/@polyadic/size.m

- Signature: `varargout=size(p,dim)`

## Purpose

Returns the represented matrix dimensions. With no `dim`, it returns `[nrows ncols]` as one output or the row and column counts as two outputs; with `dim` equal to 1 or 2, it returns that dimension.

## Physical / mathematical content

The row count is obtained from the first prefix matrix when present (otherwise from the row dimensions of the first core term); the column count is obtained from the last suffix matrix when present (otherwise from the column dimensions of the first core term).

## Numerical / algorithmic content

## Parameters / inputs

- `p`: a polyadic object
- `dim` (optional): dimension whose size is required; must be 1 or 2

## Outputs

- One output: `[nrows ncols]` or the requested dimension size; two outputs: `nrows` and `ncols`

## Implementation structure

- Computes row and column counts from the polyadic core dimensions and any prefix/suffix matrices.
- Supports `size(p)`, `[nrows,ncols]=size(p)`, and `size(p,dim)` for `dim` 1 or 2.
- Check consistency
- Get row dimension
- The leftmost matrix in the prefix
- The cores of the polyadic
- Get column dimension
- The rightmost matrix in the suffix
- Compose the answer
