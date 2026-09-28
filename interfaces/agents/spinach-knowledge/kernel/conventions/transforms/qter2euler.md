# kernel/conventions/transforms/qter2euler.m

- Signature: `[alpha,beta,gamma]=qter2euler(q)`

## Purpose

Converts a quaternion in the active convention to Euler angles in the ZYZ active convention, matching `euler2dcm.m`.

## Syntax

```matlab
[alpha,beta,gamma]=qter2euler(q)
```

## Parameters / inputs

- `q` — structure with fields `q.u`, `q.i`, `q.j`, and `q.k` containing the four quaternion components. Each field must be a real scalar or column vector; the fields must have the same number of elements. Vector inputs are converted elementwise. Quaternion norms must be significantly nonzero.

## Outputs

- `alpha`, `beta`, `gamma` — Euler angles in radians (ZYZ active convention), with the same shape as the quaternion component fields. Euler angles are not unique; the returned angles satisfy `euler2dcm(alpha,beta,gamma)=qter2dcm(q)`, with `beta` in `[0,pi]`.

## Numerical / algorithmic content

The function normalizes each quaternion before computing the angles. It calculates `sum_ag=2*atan2(q.k,q.u)` and `dif_ga=2*atan2(q.i,q.j)`, then sets `beta=2*atan2(sqrt(q.i.^2+q.j.^2),sqrt(q.u.^2+q.k.^2))`, `alpha=(sum_ag-dif_ga)/2`, and `gamma=(sum_ag+dif_ga)/2`.

[Source documentation](https://spindynamics.org/wiki/index.php?title=qter2euler.m).