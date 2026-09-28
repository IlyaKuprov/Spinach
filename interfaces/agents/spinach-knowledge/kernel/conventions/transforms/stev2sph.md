# kernel/conventions/transforms/stev2sph.m

- Signature: `Bkq=stev2sph(k,Bkq)`

## Purpose

Converts coefficients of Stevens operators produced by `stevens.m` into coefficients of irreducible spherical tensor operators produced by `irr_sph_ten.m`. Supports spherical ranks 1 through 12. Source for ranks up to 6: http://dx.doi.org/10.1088/0022-3719/18/7/009

## Syntax

```matlab
Bkq=stev2sph(k,Bkq)
```

## Parameters / inputs

- `k`: spherical rank; a real integer from 1 to 12.
- `Bkq`: a finite, real column vector of `2*k+1` coefficients in front of Stevens operators, in increasing order of projections.

## Output

- `Bkq`: a column vector of `2*k+1` complex coefficients in front of irreducible spherical tensor operators, in decreasing order of projections.

## Numerical / algorithmic content

The implementation selects a rank-specific vector of scaling factors, constructs a transformation matrix with diagonal and antidiagonal terms, and applies it to the input coefficients as `transpose(Bkq'*A)`. It checks that `k` is a real integer from 1 to 12 and that the input is a finite, real column vector of length `2*k+1`.

For ranks 7–12, the squared scaling factors are exact rationals computed from the integer coefficient table of `stevens.m` (Ryabov, J. Magn. Reson. 140, 141 (1999)) and the normalization of `irr_sph_ten.m`: `2^(k-2)*P(k,q)/C(k,q)^2` for `q>0`, and `2^k*P(k,0)/C(k,0)^2` for `q=0`. Here `P(k,q)` is the product of `(k+p)(k-p+1)` for `p` from `q+1` to `k`, and `C(k,q)` is the `stevens.m` coefficient, doubled for even `k` and odd `q`. The same expression reproduces the published ranks 1–6. The Ryabov table has a cluster of large primes at rank 9 projections 1 and 2, hence the denominators there.

## Source

- <https://spindynamics.org/wiki/index.php?title=stev2sph.m>
- e.suturina@bath.ac.uk
- ilya.kuprov@weizmann.ac.il