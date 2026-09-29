# examples/giant_spin/dy_lft_single_1.m

- MATLAB implementation: [examples/giant_spin/dy_lft_single_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/dy_lft_single_1.m)

- Source: [examples/giant_spin/dy_lft_single_1.m](../../../../../examples/giant_spin/dy_lft_single_1.m)
- Signature: `dy_lft_single_1()` (no input or output arguments)

## Purpose

Reproduces a MOLCAS ligand-field calculation for a single Dy(III) ion. The source labels the runtime as seconds; this is a source comment, not a timing measured here.

## Spin model and construction

The system contains `E16` at `sys.magnet=0`. The real anisotropic Zeeman matrix is `V' * diag(D) * V`, with principal values `D=[1.325781, 1.322640, 1.317917]`. Rank-2, rank-4, and rank-6 MOLCAS ligand-field coefficient arrays `Bkq` are converted with `icm2hz` and `stev2sph`, then transformed using Wigner matrices and the Euler angles derived from the molecular-frame rotations `R` and `rkd`. The source comment on `rkd` says its entries need more decimal places.

Spinach receives these tensors in `inter.giant.coeff`; the listed giant-spin Euler angles are zero. The basis is `zeeman-hilb` with `approximation='none'`, then the system is created and passed to `basis`.

## Calculation and output

The script evaluates `geffect(spin_system,[1 2])` and prints its eigenvalues alongside the embedded MOLCAS reference `19.2967, 0.0529, 0.0579`. These are reference values written in the source; this note does not assert a fresh run or independent agreement. The function does not return a spectrum or other output argument.

The source does not annotate units for the `Bkq` input arrays; it does explicitly pass them through `icm2hz`. No applied-field sweep is performed.
