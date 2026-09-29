# kernel/utilities/wigner.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/wigner.m)

## Purpose

Computes Wigner D matrices for a given rank `l` and three Euler angles, following the definition (Brink & Satchler, Eq. 2.13):

```
D = expm(-1i*Lz*alp) * expm(-1i*Ly*bet) * expm(-1i*Lz*gam)
```

where the generic branch takes the `(2*l+1)`-dimensional angular-momentum generators from `pauli(2*l+1)`. The ZYZ convention is used for the Euler angles (see Brink and Satchler, Figures 1 and 2).

## Behaviour

- Syntax: `D = wigner(l, alp, bet, gam)`.
- The rank `l` may be a non-negative integer or half-integer.
- For `l == 1` and `l == 2`, the matrix elements are hard-coded explicitly for speed; the source comments note these run faster than a generic `expm` call.
- For all other ranks, the function obtains Pauli matrices via `pauli(2*l+1)` and computes the matrix product of matrix exponentials directly.
- Input validation is performed by an internal `grumble` subfunction:
  - `l` must be a non-negative real scalar integer or half-integer, otherwise the function errors with `'l must be a non-negative real integer or half-integer.'`.
  - `alp`, `bet` and `gam` must be real scalars, otherwise the function errors with `'alp, bet and gam must be real scalars.'`.
- The output matrix `D` has rows and columns sorted by descending magnetic quantum numbers. For example, for `l = 2` the layout is:

```
[D( 2, 2)  ...  D( 2,-2)
   ...      ...    ...
[D(-2, 2)  ...  D(-2,-2)]
```

- The output is intended to be used as `y = D*x`, where `x` is a column vector of irreducible spherical tensor coefficients listed vertically in the order `T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2)` (for `l = 2`).

## Inputs and outputs

| Name | Type | Description |
|------|------|-------------|
| `l` | scalar | Rank of the Wigner matrix; may be half-integer. Must be a non-negative real integer or half-integer. |
| `alp` | scalar | Euler angle, radians. Must be a real scalar. |
| `bet` | scalar | Euler angle, radians. Must be a real scalar. |
| `gam` | scalar | Euler angle, radians. Must be a real scalar. |
| `D` | matrix | Wigner D matrix of size `(2l+1) x (2l+1)`, with rows and columns ordered by descending magnetic quantum number. |

## References

- Brink, D. M. and Satchler, G. R., *Angular Momentum*, Eq. 2.13 and Figures 1 and 2 (cited in the source header for the definition and the ZYZ Euler angle convention).
- Spinach Wiki page for `wigner.m`: [https://spindynamics.org/wiki/index.php?title=wigner.m](https://spindynamics.org/wiki/index.php?title=wigner.m)
