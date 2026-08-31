# kernel/conventions/transforms/anax2qter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/anax2qter.m`
- Signature: `q=anax2qter(rot_axis,rot_angle)`
- Total lines: 62

## Purpose

Converts angle-axis rotation parameters into a quaternion. Syntax: q=anax2qter(rot_axis,rot_angle)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(rot_axis,rot_angle)`.
- Lines 27-28: Normalize the axis vector; implemented by `rot_axis=rot_axis/norm(rot_axis,2)`.
- Lines 30-31: Compute the quaternion; implemented by `q.u=cos(rot_angle/2)`.

### Key state/data transformations

- Lines 28: computes `rot_axis` using `rot_axis=rot_axis/norm(rot_axis,2)`.
- Lines 31: computes `q.u` using `q.u=cos(rot_angle/2)`.
- Lines 32: computes `q.i` using `q.i=rot_axis(1)*sin(rot_angle/2)`.
- Lines 33: computes `q.j` using `q.j=rot_axis(2)*sin(rot_angle/2)`.
- Lines 34: computes `q.k` using `q.k=rot_axis(3)*sin(rot_angle/2)`.

### Local helper functions

- Line 39: `grumble()` — `function grumble(rot_axis,rot_angle)`.
  - Representative operation: `if (~isnumeric(rot_axis))||(~isnumeric(rot_angle))`.
  - Representative operation: `error('both inputs must be numeric.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rot_axis()`, `any()`.
