# kernel/conventions/transforms/spsk2mat.m

- Signature: `M=spsk2mat(iso,sp,sk,alp,bet,gam)`

## Purpose

Converts the span and skew representation of a 3x3 interaction tensor (Herzfeld–Berger convention) into a matrix. Euler angles are in radians.

## Parameters / inputs

- `iso` — isotropic part, `(xx+yy+zz)/3`, where `xx`, `yy`, and `zz` are the eigenvalues.
- `sp` — span: the largest eigenvalue minus the smallest eigenvalue.
- `sk` — skew, `3*(yy-iso)/sp`, where `yy` is the middle eigenvalue.
- `alp` — alpha Euler angle in radians.
- `bet` — beta Euler angle in radians.
- `gam` — gamma Euler angle in radians.

All inputs must be real scalars. Skew must lie in `[-1,+1]`, and span must be nonnegative.

## Numerical / algorithmic content

The eigenvalues are computed as `xx=iso-(3+sk)*sp/6`, `yy=iso+sk*sp/3`, and `zz=iso+(3-sk)*sp/6`. With `R=euler2dcm(alp,bet,gam)`, the matrix is `M=R*diag([xx yy zz])*R'`.

## Outputs

- `M` — 3x3 matrix.

The reverse transformation is ill-defined.

Source reference: <https://spindynamics.org/wiki/index.php?title=spsk2mat.m>