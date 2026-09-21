# kernel/conventions/transforms/qter2dcm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/qter2dcm.m`
- Signature: `dcm=qter2dcm(q)`
- Total lines: 60

## Purpose

Converts a unit quaternion into a direction cosine matrix in the active convention, matching the one used by euler2dcm.m function. Syntax: dcm=qter2dcm(q)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion

## Outputs

- dcm -directional cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Note: Matlab's Aerospace Toolbox quat2dcm() returns the
- transpose of this matrix for the same quaternion.

## Implementation structure

- Converts a unit quaternion into a direction cosine matrix in
- the active convention, matching the one used by euler2dcm.m
- function. Syntax:
- dcm=qter2dcm(q)
- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion
- dcm -directional cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Note: Matlab's Aerospace Toolbox quat2dcm() returns the
- transpose of this matrix for the same quaternion.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `all()`, `isfield()`, `isscalar()`, `eps()`.
