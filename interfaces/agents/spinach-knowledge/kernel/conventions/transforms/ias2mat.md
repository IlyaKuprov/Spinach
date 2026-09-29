# kernel/conventions/transforms/ias2mat.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/ias2mat.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ias2mat.m)

## Signature

`C = ias2mat(a, d, A)`

## Purpose and decomposition

Reconstructs a real 3-by-3 interaction matrix from isotropic, antisymmetric, and symmetric components. The source defines its convention for real vectors `u` and `v` by:

`a*(u'*v) + d'*cross(u,v) + u'*A*v = u'*C*v`

The constructed matrix is exactly:

`C = a*eye(3,3) + [0 d(3) -d(2); -d(3) 0 d(1); d(2) -d(1) 0] + A`

This fixes the antisymmetric-component sign/orientation; it is the matrix written in the implementation, not an inferred inverse transform.

## Inputs and output

- `a`: real numeric scalar.
- `d`: real numeric 3-by-1 column vector.
- `A`: real numeric 3-by-3 matrix. The source rejects it when `norm(A-A',2)/norm(A,2) > 1e-6`.
- `C`: reconstructed real 3-by-3 matrix.

## Reference

- [Spinach Wiki: ias2mat.m](https://spindynamics.org/wiki/index.php?title=ias2mat.m)
