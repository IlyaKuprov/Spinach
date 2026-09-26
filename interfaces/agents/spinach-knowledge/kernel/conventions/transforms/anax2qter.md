# kernel/conventions/transforms/anax2qter.m

- Signature: `q=anax2qter(rot_axis,rot_angle)`

## Purpose

Converts angle-axis rotation parameters into a quaternion. Syntax: q=anax2qter(rot_axis,rot_angle)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- rot_axis -cartesian direction vector given as a row or column
- with three real elements
- rot_angle -rotation angle in radians

## Outputs

- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion

## Implementation structure

- Converts angle-axis rotation parameters into a quaternion. Syntax:
- q=anax2qter(rot_axis,rot_angle)
- rot_axis -cartesian direction vector given as a row or column
- with three real elements
- rot_angle -rotation angle in radians
- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion
- Check consistency
- Normalize the axis vector
- Compute the quaternion
- Consistency enforcement
- Some lecture videos in IK's Spin Dynamics course (https://spindynamics.org)
