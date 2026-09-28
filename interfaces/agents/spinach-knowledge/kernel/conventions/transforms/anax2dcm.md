# kernel/conventions/transforms/anax2dcm.m

- Signature: `dcm=anax2dcm(rot_axis,rot_angle)`

## Purpose

Converts an angle and axis to a direction cosine matrix in the active convention used by `euler2dcm.m`. The angle is in radians, and the function normalizes the axis.

## Physical / mathematical content

The matrix applies an active rotation: `v=R*v` for a 3×1 vector and `A=R*A*R'` for a 3×3 interaction tensor. MATLAB Aerospace Toolbox's `quat2dcm()` returns the transpose of this matrix for the same rotation.

## Numerical / algorithmic content

The function normalizes the axis, then constructs the matrix as `eye(3) + sin(rot_angle)*K + (1-cos(rot_angle))*(rot_axis*rot_axis'-eye(3))`, where `K` is the skew-symmetric matrix formed from the axis components.

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

The function checks that both inputs are real and numeric, that `rot_axis` has three elements and is nonzero, and that `rot_angle` is scalar. It then normalizes the axis and computes the direction cosine matrix.
