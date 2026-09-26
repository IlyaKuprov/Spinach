# kernel/conventions/transforms/qter2anax.m

- Signature: `[rot_axis,rot_angle]=qter2anax(q)`

## Purpose

Converts a quaternion representation of a rotation into angle-axis rotation parameters. Syntax: [rot_axis,rot_angle]=qter2anax(q)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- q -quaternion, a structure with four fields
- q.u, q.i, q.j, q.k giving the four compo-
- nents of the quaternion

## Outputs

- rot_axis -cartesian direction vector as a row
- with three real elements
- rot_angle -rotation angle in radians

## Implementation structure

- Converts a quaternion representation of a rotation into angle-axis
- rotation parameters. Syntax:
- [rot_axis,rot_angle]=qter2anax(q)
- q - quaternion, a structure with four fields
- q.u, q.i, q.j, q.k giving the four compo-
- nents of the quaternion
- rot_axis -cartesian direction vector as a row
- with three real elements
- rot_angle -rotation angle in radians
- Check consistency
- Normalize the quaternion
- Compute the vector part norm
