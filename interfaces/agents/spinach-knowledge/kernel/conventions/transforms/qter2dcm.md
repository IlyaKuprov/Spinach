# kernel/conventions/transforms/qter2dcm.m

- Signature: `dcm=qter2dcm(q)`

## Purpose

Converts a unit quaternion into a direction cosine matrix in the active convention, matching the one used by euler2dcm.m function. Syntax: dcm=qter2dcm(q)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion

## Outputs

- dcm -directional cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Note: Matlab's Aerospace Toolbox quat2dcm() returns the
- transpose of this matrix for the same quaternion.

## Implementation structure

- Converts a unit quaternion into a direction cosine matrix in
- the active convention, matching the one used by euler2dcm.m
- function. Syntax:
- dcm=qter2dcm(q)
- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion
- dcm -directional cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Note: Matlab's Aerospace Toolbox quat2dcm() returns the
- transpose of this matrix for the same quaternion.
