# kernel/conventions/transforms/qter2anax.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/qter2anax.m`
- Signature: `[rot_axis,rot_angle]=qter2anax(q)`
- Total lines: 60

## Purpose

Converts a quaternion representation of a rotation into angle-axis rotation parameters. Syntax: [rot_axis,rot_angle]=qter2anax(q)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `atan2()`, `all()`, `isfield()`.
