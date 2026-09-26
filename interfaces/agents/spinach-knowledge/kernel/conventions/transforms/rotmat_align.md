# kernel/conventions/transforms/rotmat_align.m

- Signature: `rot_mat=rotmat_align(v_from,v_to)`

## Purpose

Rotation matrix aligning one vector with another vector. Syntax: rot_mat=rotmat_align(v_from,v_to)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- v_from -three-element real vector to rotate
- v_to -three-element real vector to align to

## Outputs

- rot_mat -3x3 rotation matrix that satisfies
- rot_mat*(v_from/norm(v_from,2))=v_to/norm(v_to,2)
- Note: aligning one vector with another leaves one rotational degree of
- freedom around the aligned direction. If the resulting matrix is
- converted to ZYZ Euler angles, this freedom appears as a non-uni-
- que third Euler angle (the twist around the aligned axis). This
- implementation fixes that freedom by returning the minimum-
- angle alignment without any additional twist around the aligned
- direction. In the anti-parallel case the axis is non-unique, and
- the first null-space basis vector orthogonal to v_from is used.

## Implementation structure

- Rotation matrix aligning one vector with another vector. Syntax:
- rot_mat=rotmat_align(v_from,v_to)
- v_from -three-element real vector to rotate
- v_to -three-element real vector to align to
- rot_mat -3x3 rotation matrix that satisfies
- rot_mat*(v_from/norm(v_from,2))=v_to/norm(v_to,2)
- Note: aligning one vector with another leaves one rotational degree of
- freedom around the aligned direction. If the resulting matrix is
- converted to ZYZ Euler angles, this freedom appears as a non-uni-
- que third Euler angle (the twist around the aligned axis). This
- implementation fixes that freedom by returning the minimum-
- angle alignment without any additional twist around the aligned
