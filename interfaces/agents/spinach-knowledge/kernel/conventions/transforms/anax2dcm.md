# kernel/conventions/transforms/anax2dcm.m

- Signature: `dcm=anax2dcm(rot_axis,rot_angle)`

## Purpose

Converts angle-axis rotation parameters to a direction cosine matrix in the active convention, matching the one used by euler2dcm.m function. Angle should be in radians, axis is normalized by the function. Syntax: dcm=anax2dcm(rot_axis,rot_angle)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- rot_axis -cartesian direction vector given as
- a row or column with three real ele-
- ments
- rot_angle -rotation angle in radians

## Outputs

- dcm -directional cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Note: Matlab's Aerospace Toolbox quat2dcm() returns the
- transpose of this matrix for the same rotation.

## Implementation structure

- Converts angle-axis rotation parameters to a direction
- cosine matrix in the active convention, matching the one
- used by euler2dcm.m function. Angle should be in radians,
- axis is normalized by the function. Syntax:
- dcm=anax2dcm(rot_axis,rot_angle)
- rot_axis -cartesian direction vector given as
- a row or column with three real ele-
- ments
- rot_angle -rotation angle in radians
- dcm -directional cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
