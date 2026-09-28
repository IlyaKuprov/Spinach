# kernel/conventions/transforms/qform2sph.m

- Signature: `[r0,r1,r2]=qform2sph(A)`

## Purpose

Expands the normalized quadratic form `[x y z]*A*[x y z]'/ (x^2+y^2+z^2)` in spherical harmonics, returning coefficients for ranks 0, 1, and 2.

## Physical / mathematical content

The rank-0 coefficient is `r0=(2/3)*sqrt(pi)*trace(A)`. Rank-1 coefficients are zero. The five rank-2 coefficients are ordered by m=2,1,0,-1,-2; the source evaluates them from the symmetric matrix elements.

## Numerical / algorithmic content

The function requires a real numeric symmetric 3x3 matrix. It returns a scalar rank-0 coefficient, three zero rank-1 coefficients, and five rank-2 coefficients. For the rank-2 terms, let `c=sqrt(2*pi/15)`; their values in order are `c*(A11-A22-2i*A12)`, `-2*c*(A13-i*A23)`, `(2/3)*sqrt(pi/5)*(2*A33-A11-A22)`, `2*c*(A13+i*A23)`, and `c*(A11-A22+2i*A12)`.

## Parameters / inputs

- `A` — real numeric symmetric 3x3 matrix.

## Outputs

- `r0` — rank-0 coefficient.
- `r1` — three zero rank-1 coefficients.
- `r2` — five rank-2 coefficients, ordered by m=2,1,0,-1,-2.

## Implementation structure

The function checks that A is real, numeric, symmetric, and 3x3, then evaluates the rank-0 and rank-2 expressions; rank 1 is set to zero.
