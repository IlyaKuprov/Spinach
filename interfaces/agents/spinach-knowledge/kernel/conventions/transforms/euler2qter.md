# kernel/conventions/transforms/euler2qter.m

- Signature: `q=euler2qter(arg1,arg2,arg3)`

## Purpose

Converts Euler angles (ZYZ active convention) into a unit quaternion in the active convention, matching euler2dcm.m function. Syntax: q=euler2qter(alpha,beta,gamma) OR q=euler2qter([alpha beta gamma])

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), scalars or column vec-
- tors of equal length

## Outputs

- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; for column
- vector inputs each field is a column vector
- Note: the quaternion returned represents the same rotation
- as euler2dcm(alpha,beta,gamma); it is converted into
- that matrix by qter2dcm.m function.

## Implementation structure

- Converts Euler angles (ZYZ active convention) into a unit
- quaternion in the active convention, matching euler2dcm.m
- function. Syntax:
- q=euler2qter(alpha,beta,gamma)
- q=euler2qter([alpha beta gamma])
- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), scalars or column vec-
- tors of equal length
- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; for column
- vector inputs each field is a column vector
- Note: the quaternion returned represents the same rotation
