# kernel/utilities/wigner.m

- Signature: `D=wigner(l,alp,bet,gam)`

## Purpose

Computes a Wigner D matrix using ZYZ Euler angles (Brink and Satchler, Eq. 2.13; Figures 1 and 2). The matrix represents `expm(-1i*L.z*alp)*expm(-1i*L.y*bet)*expm(-1i*L.z*gam)`, where `L` is obtained from `pauli(2*l+1)`.

## Parameters / inputs

- `l`: non-negative integer or half-integer rank.
- `alp`, `bet`, `gam`: real scalar Euler angles in radians.

## Output

- `D`: Wigner D matrix with rows and columns ordered by descending magnetic quantum number, from `l` to `-l`. For `l=2`, the first row runs from `D(2,2)` to `D(2,-2)`, and the last from `D(-2,2)` to `D(-2,-2)`. Apply it as `y=D*x` to a column of irreducible spherical tensor coefficients ordered `T(2,2)`, `T(2,1)`, `T(2,0)`, `T(2,-1)`, `T(2,-2)`.

## Implementation

Inputs are checked for the stated types and ranges. Ranks `l=1` and `l=2` use hard-coded matrices for speed; other ranks use the product of three matrix exponentials above.

Source: <https://spindynamics.org/wiki/index.php?title=wigner.m>

Contact: ilya.kuprov@weizmann.ac.il