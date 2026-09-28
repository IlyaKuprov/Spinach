# kernel/conventions/transforms/tsm2param.m

- Signature: `[ax,rh,angles]=tsm2param(M)`

## Purpose

Converts a traceless symmetric 3x3 interaction matrix to axiality, rhombicity, and Euler angles. The source warns that this conversion is unstable and recommends publishing the 3x3 matrix instead, citing IUPAC.

## Parameters / inputs

- `M` — a 3x3 matrix or five independent elements ordered `[Mxx, Mxy, Mxz, Myy, Myz]`. For five elements, the matrix is assembled symmetrically with `Mzz = -Mxx - Myy`.

## Outputs

- `ax` — axiality: `2*DZ - (DX + DY)`.
- `rh` — rhombicity: `DY - DX`.
- `angles` — Euler angles in radians, one of eight equivalent sets.

Eigenvalues use Mehring ordering: `DZ` is the largest and `DX` the smallest, including signs; `DY` is the remaining eigenvalue. The eigenvectors are ordered X, Y, Z, adjusted to have positive determinant, and passed to `dcm2euler`. The source defines a consistency-checking helper, but the call to it is commented out.

[Source](https://spindynamics.org/wiki/index.php?title=tsm2param.m)