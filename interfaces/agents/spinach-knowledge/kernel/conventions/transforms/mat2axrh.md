# kernel/conventions/transforms/mat2axrh.m

MATLAB source: [kernel/conventions/transforms/mat2axrh.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/mat2axrh.m)
Spinach Wiki: [mat2axrh.m](https://spindynamics.org/wiki/index.php?title=mat2axrh.m)

## Purpose and usage

Extracts isotropic value, axiality, rhombicity, and eigenvalues from a real symmetric 3-by-3 interaction matrix.

`[iso, ax, rh, eigvals] = mat2axrh(M)`

## Input and outputs

- `M`: a real, symmetric 3-by-3 matrix.
- `iso`: scalar isotropic component, the mean of the three eigenvalues.
- `ax`: scalar axiality.
- `rh`: scalar rhombicity.
- `eigvals`: the three eigenvalues sorted in ascending Mehring order, as `[xx; yy; zz]` with `xx <= yy <= zz`.

The implementation checks that `M` is numeric, real, 3-by-3, and symmetric. It does not return Euler angles: the corresponding transformation is ill-defined.

## Tensor convention

If the eigenvalues in Mehring order are `xx <= yy <= zz`, the current definitions are

- `iso = (xx + yy + zz) / 3`;
- `ax = 2*zz - (xx + yy)`;
- `rh = yy - xx`.

These are the present definitions used by `axrh2mat.m`; calling that function with `iso`, `ax`, and `rh` reproduces the eigenvalues of `M`. Earlier versions of `mat2axrh` used `ax = zz - (xx + yy)/2` and `rh = xx - yy`. Those historical definitions differ by a factor of two in axiality and a sign reversal in rhombicity; they are not the current outputs.
