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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(q)`.
- Lines 33-34: Normalise the quaternion; implemented by `qnorm=norm([q.u q.i q.j q.k],2)`.
- Lines 38-39: Compute the DCM; implemented by `dcm=[1-2*(q.j^2+q.k^2) 2*(q.i*q.j-q.u*q.k) 2*(q.i*q.k+q.u*q.j)`.

### Key state/data transformations

- Lines 34: computes `qnorm` using `qnorm=norm([q.u q.i q.j q.k],2)`.
- Lines 35: computes `q.u` using `q.u=q.u/qnorm; q.i=q.i/qnorm`.
- Lines 36: computes `q.j` using `q.j=q.j/qnorm; q.k=q.k/qnorm`.
- Lines 39: computes `dcm` using `dcm=[1-2*(q.j^2+q.k^2) 2*(q.i*q.j-q.u*q.k) 2*(q.i*q.k+q.u*q.j)`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(q)`.
  - Representative operation: `if ~all(isfield(q,{'i','j','k','u'}))`.
  - Representative operation: `error('quaternion data structure must contain u, i, j, and k fields.')`.

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
