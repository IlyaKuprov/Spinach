# kernel/conventions/transforms/anax2dcm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/anax2dcm.m`
- Signature: `dcm=anax2dcm(rot_axis,rot_angle)`
- Total lines: 70

## Purpose

Converts angle-axis rotation parameters to a direction cosine matrix in the active convention, matching the one used by euler2dcm.m function. Angle should be in radians, axis is normalized by the function. Syntax: dcm=anax2dcm(rot_axis,rot_angle)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(rot_axis,rot_angle)`.
- Lines 37-38: Normalize the axis; implemented by `rot_axis=rot_axis(:)/norm(rot_axis(:),2)`.
- Lines 40-41: Compute the DCM; implemented by `dcm=eye(3)+sin(rot_angle)*[ 0 -rot_axis(3) rot_axis(2)`.

### Key state/data transformations

- Lines 38: computes `rot_axis` using `rot_axis=rot_axis(:)/norm(rot_axis(:),2)`.
- Lines 41: computes `dcm` using `dcm=eye(3)+sin(rot_angle)*[ 0 -rot_axis(3) rot_axis(2)`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(rot_axis,rot_angle)`.
  - Representative operation: `if (~isnumeric(rot_axis))||(~isnumeric(rot_angle))`.
  - Representative operation: `error('both inputs must be numeric.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rot_axis()`, `any()`.
