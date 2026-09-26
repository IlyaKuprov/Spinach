# kernel/conventions/transforms/euler_sup.m

- Signature: `rot_cmp=euler_sup(rot_one,rot_two)`

## Purpose

Superposition of ZYZ active Euler rotations. Syntax: rot_cmp=euler_sup(rot_one,rot_two)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

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
