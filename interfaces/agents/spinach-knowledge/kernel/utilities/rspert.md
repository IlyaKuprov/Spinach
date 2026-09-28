# kernel/utilities/rspert.m

- Signature: `[Ep,Vp]=rspert(E0,H1,order)`

## Purpose

Computes Rayleigh-Schrodinger perturbation theory to the specified order, following Eqs. 2.21–2.23 of Stefan Stoll’s PhD thesis, with the typo in the numerator of Eq. 2.21 corrected.

## Parameters / inputs

- `E0` — Real column vector of eigenvalues of `H0`.
- `H1` — Hermitian perturbation written in the eigenbasis of `H0`.
- `order` — Positive integer specifying the perturbation order. Sixth order is the sensible maximum.

## Outputs

- `Ep` — Real vector of eigenvalues of `H0+H1` to the specified order. Its entries are not necessarily sorted in the same way as the input.
- `Vp` — Normalized eigenvectors of `H0+H1` to the specified order, returned as a square unitary matrix. Its columns correspond, in order, to the eigenvalues in `Ep`.

## Notes

`H0` must have no degenerate energy levels. The perturbation theory converges only when `norm(H1,2)` is much smaller than the smallest energy gap in `H0`. Numerical artifacts appear beyond sixth order. Computational complexity is linear in `order` and cubic in the matrix dimension.

Source: [rspert.m](https://spindynamics.org/wiki/index.php?title=rspert.m). Contact: ilya.kuprov@weizmann.ac.il.