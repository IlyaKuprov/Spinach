# kernel/conventions/transforms/qter2euler.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/qter2euler.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=qter2euler.m)

- Signature: `[alpha,beta,gamma]=qter2euler(q)`

## Purpose

Convert quaternion components in the active convention to ZYZ active Euler angles, in the convention matched to `euler2dcm.m`. The source header describes a unit quaternion; the implementation normalises each accepted quaternion before evaluating the angles.

## Inputs and checks

- `q` is a structure with numeric, real, column fields `u`, `i`, `j`, and `k`. Each field may be a scalar or column vector; all four must have the same number of elements.
- Each quaternion norm must be at least `sqrt(eps())`; smaller norms produce an error.

## Conversion

For each element, let `qhat` be `(u,i,j,k)/sqrt(u^2+i^2+j^2+k^2)`. The implementation then computes:

```text
sum_ag = 2*atan2(qhat.k,qhat.u)
dif_ga = 2*atan2(qhat.i,qhat.j)
beta   = 2*atan2(sqrt(qhat.i^2+qhat.j^2),sqrt(qhat.u^2+qhat.k^2))
alpha  = (sum_ag-dif_ga)/2
gamma  = (sum_ag+dif_ga)/2
```

`alpha`, `beta`, and `gamma` are in radians; `beta` is in `[0,pi]`. Outputs have the same scalar or column-vector shape as the quaternion component fields. Euler angles are not unique; the source documents that the returned angles satisfy `euler2dcm(alpha,beta,gamma)=qter2dcm(q)`.
