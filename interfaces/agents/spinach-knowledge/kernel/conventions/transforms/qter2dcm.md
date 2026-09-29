# kernel/conventions/transforms/qter2dcm.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/qter2dcm.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=qter2dcm.m)

- Signature: `dcm=qter2dcm(q)`

## Behaviour

The function normalises quaternion components `(u,i,j,k)` and returns the active-convention direction cosine matrix used by `euler2dcm.m`:

`[1-2*(j^2+k^2), 2*(i*j-u*k), 2*(i*k+u*j); 2*(i*j+u*k), 1-2*(i^2+k^2), 2*(j*k-u*i); 2*(i*k-u*j), 2*(j*k+u*i), 1-2*(i^2+j^2)]`

Use it on a column vector as `v=dcm*v`, and on a 3x3 interaction tensor as `A=dcm*A*dcm'`. For the same quaternion, MATLAB Aerospace Toolbox `quat2dcm()` returns the transpose of this matrix.

## Input and output

The input structure must contain numeric, real scalar fields `u`, `i`, `j`, and `k`. A Euclidean quaternion norm below `sqrt(eps())` raises an error; otherwise the components are normalised before constructing the matrix. The output `dcm` is a 3x3 direction cosine matrix.
