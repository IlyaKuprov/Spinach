# kernel/conventions/transforms/euler_equiv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/euler_equiv.m`
- Signature: `answer=euler_equiv(eulers_a,eulers_b,tol)`
- Total lines: 79

## Purpose

Checks whether two ZYZ active Euler angle sets specify the same rotation. Syntax: answer=euler_equiv(eulers_a,eulers_b,tol)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(eulers_a,eulers_b,tol)`.
- Lines 34-35: Cast Euler angles to row vectors; implemented by `eulers_a=eulers_a(:)'`.
- Lines 38-39: Build direction cosine matrices; implemented by `dcm_a=euler2dcm(eulers_a)`.
- Lines 42-43: Compute the relative rotation; implemented by `dcm_rel=dcm_b*dcm_a'`.
- Lines 45-48: Compute the sine and cosine of the relative rotation angle; implemented by `sin_ang=norm([dcm_rel(3,2)-dcm_rel(2,3); dcm_rel(1,3)-dcm_rel(3,1); dcm_rel(2,1)-dcm_rel(1,2)],2)/2`.
- Lines 51-52: Clamp round-off in the cosine; implemented by `cos_ang=max(min(cos_ang,1),-1)`.
- Lines 54-55: Compare the geodesic distance on SO(3); implemented by `answer=(atan2(sin_ang,cos_ang)<=tol)`.

### Key state/data transformations

- Lines 35: computes `eulers_a` using `eulers_a=eulers_a(:)'`.
- Lines 36: computes `eulers_b` using `eulers_b=eulers_b(:)'`.
- Lines 39: computes `dcm_a` using `dcm_a=euler2dcm(eulers_a)`.
- Lines 40: computes `dcm_b` using `dcm_b=euler2dcm(eulers_b)`.
- Lines 43: computes `dcm_rel` using `dcm_rel=dcm_b*dcm_a'`.
- Lines 46-48: computes `sin_ang` using `sin_ang=norm([dcm_rel(3,2)-dcm_rel(2,3); dcm_rel(1,3)-dcm_rel(3,1); dcm_rel(2,1)-dcm_rel(1,2)],2)/2`.
- Lines 49: computes `cos_ang` using `cos_ang=(trace(dcm_rel)-1)/2`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(eulers_a,eulers_b,tol)`.
  - Representative operation: `if (~isnumeric(eulers_a))||(~isreal(eulers_a))||(~isvector(eulers_a))|| (numel(eulers_a)~=3)||any(~isfinite(eulers_a(:)))`.
  - Representative operation: `(numel(eulers_a)~=3)||any(~isfinite(eulers_a(:)))`.

## Parameters / inputs

- eulers_a -first Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- eulers_b -second Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- tol -non-negative angular tolerance, radians

## Outputs

- answer -true if the relative rotation angle between
- the two rotations is not greater than tol
- Note: Euler angles are not unique, and so this function compares
- the rotations produced by euler2dcm(), not the angles
- themselves.

## Implementation structure

- Checks whether two ZYZ active Euler angle sets specify the same
- rotation. Syntax:
- answer=euler_equiv(eulers_a,eulers_b,tol)
- eulers_a -first Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- eulers_b -second Euler angle set [alpha beta gamma],
- tol -non-negative angular tolerance, radians
- answer -true if the relative rotation angle between
- the two rotations is not greater than tol
- Note: Euler angles are not unique, and so this function compares
- the rotations produced by euler2dcm(), not the angles
- themselves.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `eulers_a()`, `eulers_b()`, `euler2dcm()`, `dcm_rel()`, `atan2()`, `isvector()`, `any()`.
