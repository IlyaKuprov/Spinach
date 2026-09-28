# kernel/conventions/transforms/rotmat_align.m

- Signature: `rot_mat=rotmat_align(v_from,v_to)`

## Purpose

Return a 3x3 rotation matrix that aligns `v_from` with `v_to`.

## Parameters / inputs

- `v_from`: three-element real, finite vector to rotate.
- `v_to`: three-element real, finite vector to align to.
- Neither vector may have a 2-norm below `eps('double')`.

## Outputs

- `rot_mat`: 3x3 rotation matrix satisfying `rot_mat*(v_from/norm(v_from,2))=v_to/norm(v_to,2)`.

## Numerical / algorithmic content

The function normalizes both inputs, then uses their cross product as the rotation axis and their dot product as the cosine of the alignment angle. For non-collinear vectors, it constructs the matrix with Rodrigues' rotation formula. When the axis norm is below `1e-12`, it returns the identity for parallel vectors; for anti-parallel vectors, it uses the first null-space basis vector orthogonal to `v_from` as the rotation axis.

Aligning two vectors leaves a rotational degree of freedom around the aligned direction. The implementation uses the minimum-angle alignment without an additional twist; the anti-parallel rotation axis is non-unique.

Source reference: <https://spindynamics.org/wiki/index.php?title=rotmat_align.m>