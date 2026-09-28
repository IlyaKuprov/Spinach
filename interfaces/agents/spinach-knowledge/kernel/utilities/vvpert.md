# kernel/utilities/vvpert.m

- Signature: `[Ep,G]=vvpert(E0,H1,order)`

## Purpose

Computes Van Vleck perturbation theory following Shavitt and Redmon, excluding the quasi-degenerate split.

## Parameters

- `E0` — eigenvalues of `H0`, as a real column vector.
- `H1` — perturbation in the basis that diagonalises `H0`; must be finite and Hermitian.
- `order` — positive integer specifying the perturbation order. Numerical artefacts typically appear beyond order 10–12.

## Outputs

- `Ep` — real column vector of eigenvalues of `H0+H1` to the specified order; not necessarily sorted in the same way as `E0`.
- `G` — Van Vleck transformation generator. `expm(G)` is a square unitary matrix whose columns are eigenvectors in the order of `Ep`.

## Notes

`H0` must have no degenerate energy levels. The theory converges only when `norm(H1,2)` is much smaller than the smallest energy gap in `H0`. Computational complexity is cubic in both `order` and matrix dimension.

Source: <https://spindynamics.org/wiki/index.php?title=vvpert.m>

Contact: ilya.kuprov@weizmann.ac.il