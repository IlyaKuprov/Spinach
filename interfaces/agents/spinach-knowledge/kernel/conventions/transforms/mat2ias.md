# kernel/conventions/transforms/mat2ias.m

MATLAB source: [kernel/conventions/transforms/mat2ias.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/mat2ias.m)
Spinach Wiki: [mat2ias.m](https://spindynamics.org/wiki/index.php?title=mat2ias.m)

## Purpose and usage

Decomposes a real 3-by-3 interaction matrix `C` into isotropic, antisymmetric, and traceless symmetric parts:

`[a, d, A] = mat2ias(C)`

For real vectors `u` and `v`, the decomposition obeys

`u' * C * v = a * (u' * v) + d' * cross(u, v) + u' * A * v`.

## Input and outputs

- `C`: any real numeric 3-by-3 matrix; symmetry is not required.
- `a`: scalar isotropic component, `trace(C)/3`.
- `d`: 3-by-1 antisymmetric coupling vector, `[(C(2,3)-C(3,2)); (C(3,1)-C(1,3)); (C(1,2)-C(2,1))]/2`.
- `A`: 3-by-3 symmetric traceless matrix, `(C + C')/2 - a*eye(3,3)`.

The implementation checks that `C` is numeric, real, and exactly 3-by-3. It does not require `C` to be symmetric.
