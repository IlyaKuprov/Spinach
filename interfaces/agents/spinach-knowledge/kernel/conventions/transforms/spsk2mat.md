# kernel/conventions/transforms/spsk2mat.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/spsk2mat.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=spsk2mat.m)

- Signature: `M=spsk2mat(iso,sp,sk,alp,bet,gam)`

## Purpose and parameters

Convert span/skew parameters in the Herzfeld–Berger convention and three Euler angles into a `3x3` interaction matrix. The source defines `iso=(xx+yy+zz)/3`, `sp` as largest minus smallest eigenvalue, and `sk=3*(yy-iso)/sp` where `yy` is the middle eigenvalue. `alp`, `bet`, and `gam` are ZYZ active Euler angles in radians, as used by `euler2dcm.m`.

## Inputs and checks

All six inputs must be numeric, real scalars. The source rejects `abs(sk)>1` and `sp<0`; thus its explicit span guard allows zero. It does not add an explicit finite-value check.

## Eigenvalues and rotation

The principal values are calculated as:

```text
xx = iso-(3+sk)*sp/6
yy = iso+sk*sp/3
zz = iso+(3-sk)*sp/6
```

With `R=euler2dcm(alp,bet,gam)`, where `euler2dcm` uses the ZYZ active convention, the source returns `M=R*diag([xx yy zz])*R'`. Output shape is `3x3`. The source documentation notes that the reverse transformation is ill-defined; no inverse formula is supplied here.
