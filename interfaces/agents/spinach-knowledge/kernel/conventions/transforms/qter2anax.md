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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(q)`.
- Lines 28-29: Normalize the quaternion; implemented by `qnorm=norm([q.u q.i q.j q.k],2)`.
- Lines 33-34: Compute the vector part norm; implemented by `vec_norm=norm([q.i q.j q.k],2)`.
- Lines 36-37: Compute the angle and the axis; implemented by `if vec_norm==0`.

### Control flow inferred from the code

- Line 37: conditional branch on `vec_norm==0`.

### Key state/data transformations

- Lines 29: computes `qnorm` using `qnorm=norm([q.u q.i q.j q.k],2)`.
- Lines 30: computes `q.u` using `q.u=q.u/qnorm; q.i=q.i/qnorm`.
- Lines 31: computes `q.j` using `q.j=q.j/qnorm; q.k=q.k/qnorm`.
- Lines 34: computes `vec_norm` using `vec_norm=norm([q.i q.j q.k],2)`.
- Lines 38: computes `rot_angle` using `rot_angle=0; rot_axis=[0 0 1]`.
- Lines 41: computes `rot_axis` using `rot_axis=[q.i q.j q.k]/vec_norm`.

### Local helper functions

- Line 47: `grumble()` — `function grumble(q)`. "Mediocrity" does not mean an average intelligence; it means an average intelligence that resents and envies its betters.
  - Representative operation: `if ~all(isfield(q,{'i','j','k','u'}))`.
  - Representative operation: `error('quaternion data structure must contain u, i, j, and k fields.')`.

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
