# kernel/conventions/transforms/zfs2mat.m

- Signature: `M=zfs2mat(D,E,alp,bet,gam)`

## Purpose

Converts the zero-field splitting parameters `D` and `E` into a spin interaction matrix, following the convention described in the abstract of [doi:10.1063/1.1682294](http://dx.doi.org/10.1063/1.1682294).

## Physical / mathematical content

The function constructs the zero-field splitting tensor in its eigenframe as a diagonal matrix with entries `-D/3+E`, `-D/3-E`, and `2*D/3`. It then rotates the tensor using the direction-cosine matrix returned by `euler2dcm(alp,bet,gam)`. The output is made symmetric and traceless to remove floating-point residuals.

## Numerical / algorithmic content

The implementation checks that all five inputs are real numeric scalars, builds the diagonal tensor, applies `M=R*M*R'` with `R=euler2dcm(alp,bet,gam)`, and symmetrizes and removes the trace from the result. Angles are passed to `euler2dcm` in radians.

## Parameters / inputs

- `D`, `E` — real scalar zero-field splitting parameters, in Hz.
- `alp`, `bet`, `gam` — alpha, beta, and gamma Euler angles, in radians.

## Outputs

- `M` — symmetric 3×3 spin interaction matrix, in Hz.

## References

- [Zero-field splitting convention](http://dx.doi.org/10.1063/1.1682294).
- [Spin Dynamics Wiki: zfs2mat.m](https://spindynamics.org/wiki/index.php?title=zfs2mat.m)
