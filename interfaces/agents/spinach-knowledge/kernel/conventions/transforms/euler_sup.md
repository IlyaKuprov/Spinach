# kernel/conventions/transforms/euler_sup.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/euler_sup.m`
- Signature: `rot_cmp=euler_sup(rot_one,rot_two)`
- Total lines: 107

## Purpose

Superposition of ZYZ active Euler rotations. Syntax: rot_cmp=euler_sup(rot_one,rot_two)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- rot_one -first Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- rot_two -second Euler angle set [alpha beta gamma],
- radians, ZYZ active convention

## Outputs

- rot_cmp -row vector [alpha beta gamma] of the
- composite rotation, radians
- Note: rotations are applied in the supplied order
- v_rot=R_two*R_one*v
- therefore the composite matrix is
- R_comp=R_two*R_one

## Implementation structure

- Superposition of ZYZ active Euler rotations. Syntax:
- rot_cmp=euler_sup(rot_one,rot_two)
- rot_one -first Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- rot_two -second Euler angle set [alpha beta gamma],
- rot_cmp -row vector [alpha beta gamma] of the
- composite rotation, radians
- Note: rotations are applied in the supplied order
- v_rot=R_two*R_one*v
- therefore the composite matrix is
- R_comp=R_two*R_one
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rot_one()`, `rot_two()`, `wrapToPi()`, `euler2dcm()`, `all()`, `elseif()`, `dcm_comp()`, `acos()`, `dcm2euler()`, `isvector()`, `any()`.
