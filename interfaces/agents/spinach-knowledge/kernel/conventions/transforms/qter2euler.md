# kernel/conventions/transforms/qter2euler.m

- Signature: `[alpha,beta,gamma]=qter2euler(q)`

## Purpose

Converts a unit quaternion in the active convention into Euler angles (ZYZ active convention), matching euler2dcm.m function.

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Syntax

```matlab
[alpha,beta,gamma]=qter2euler(q)
```

## Parameters / inputs

- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; each field may
- be a column vector, in which case the conversion is
- performed elementwise

## Outputs

- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), same shape as the qua-
- ternion component fields
- Note: Euler angles are not unique; the angles returned sa-
- tisfy euler2dcm(alpha,beta,gamma)=qter2dcm(q) with
- beta in the [0,pi] interval.

## Implementation structure

- Converts a unit quaternion in the active convention into Euler
- angles (ZYZ active convention), matching euler2dcm.m function.
- [alpha,beta,gamma]=qter2euler(q)
- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; each field may
- be a column vector, in which case the conversion is
- performed elementwise
- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), same shape as the qua-
- ternion component fields
- Note: Euler angles are not unique; the angles returned sa-
- tisfy euler2dcm(alpha,beta,gamma)=qter2dcm(q) with
