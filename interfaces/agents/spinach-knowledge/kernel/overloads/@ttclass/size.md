# kernel/overloads/@ttclass/size.m

- Signature: `varargout=size(tt,dim)`

## Purpose

Return the row and column dimensions of the matrix represented by a tensor train, following the supported forms of MATLAB's `size` function.

## Syntax

- `sz=size(tt)` returns `[m n]`.
- `[m,n]=size(tt)` returns the row and column dimensions separately.
- `d=size(tt,dim)` returns the row dimension for `dim=1` or the column dimension for `dim=2`.

## Parameters / inputs

- `tt` — tensor-train representation of a matrix.
- `dim` — optional dimension selector, 1 or 2.

## Outputs

- `m,n` — integer dimensions of the represented matrix.

## Implementation

The function multiplies the second and third physical dimensions, respectively, across all cores. It errors for unsupported call syntax and if either resulting dimension exceeds MATLAB's `intmax`.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/size.m`](https://spindynamics.org/wiki/index.php?title=ttclass/size.m).
