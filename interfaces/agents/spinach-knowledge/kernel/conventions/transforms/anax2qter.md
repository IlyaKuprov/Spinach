# kernel/conventions/transforms/anax2qter.m

- Signature: `q=anax2qter(rot_axis,rot_angle)`

## Purpose

Converts angle-axis rotation parameters into a quaternion.

## Physical / mathematical content

The rotation axis is normalized before the quaternion is computed. For the normalized axis `a` and rotation angle `theta`, the components are `q.u=cos(theta/2)` and `[q.i,q.j,q.k]=a*sin(theta/2)`.

## Numerical / algorithmic content

The function checks that both inputs are numeric and real, that the axis has three elements and is non-zero, and that the angle is a scalar. It then normalizes the axis using its 2-norm and computes the four quaternion components.

## Parameters / inputs

- rot_axis -cartesian direction vector given as a row or column
- with three real elements
- rot_angle -rotation angle in radians

## Outputs

- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion

## Implementation structure

- `grumble` enforces the input requirements before normalization and quaternion construction.
- The source includes a comment about lecture videos in IK's Spin Dynamics course (https://spindynamics.org).
