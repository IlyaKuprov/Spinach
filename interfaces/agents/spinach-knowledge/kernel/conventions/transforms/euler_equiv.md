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
