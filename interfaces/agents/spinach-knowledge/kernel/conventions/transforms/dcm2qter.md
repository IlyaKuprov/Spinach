# kernel/conventions/transforms/dcm2qter.m

- Signature: `q=dcm2qter(dcm)`

## Purpose

Converts a direction cosine matrix in the active convention of euler2dcm.m function into a unit quaternion. Syntax: q=dcm2qter(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- dcm -directional cosine matrix, a 3x3 orthogonal
- matrix with unit determinant

## Outputs

- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion,
- normalised to q.u greater than or equal to zero
- Note: quaternions double-cover rotations; of the two candi-
- dates q and -q this function returns the one with the
- non-negative scalar part.

## Implementation structure

- Converts a direction cosine matrix in the active convention of
- euler2dcm.m function into a unit quaternion. Syntax:
- q=dcm2qter(dcm)
- dcm -directional cosine matrix, a 3x3 orthogonal
- matrix with unit determinant
- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion,
- normalised to q.u greater than or equal to zero
- Note: quaternions double-cover rotations; of the two candi-
- dates q and -q this function returns the one with the
- non-negative scalar part.
- Check consistency
