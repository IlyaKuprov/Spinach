# kernel/conventions/transforms/rotmat_align.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/rotmat_align.m`
- Signature: `rot_mat=rotmat_align(v_from,v_to)`
- Total lines: 107

## Purpose

Rotation matrix aligning one vector with another vector. Syntax: rot_mat=rotmat_align(v_from,v_to)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(v_from,v_to)`.
- Lines 34-35: Normalise the input vectors; implemented by `v_from=v_from(:)/norm(v_from,2)`.
- Lines 38-39: Set the collinearity tolerance; implemented by `align_tol=1e-12`.
- Lines 41-42: Compute the alignment axis and the alignment angle cosine; implemented by `rot_axis=cross(v_from,v_to)`.
- Lines 46-47: Handle parallel and anti-parallel vectors; implemented by `if axis_norm<align_tol`.
- Lines 50-51: Return identity when no rotation is needed; implemented by `rot_mat=eye(3)`.
- Lines 56-57: Choose an arbitrary axis orthogonal to the input vector; implemented by `null_basis=null(v_from.')`.
- Lines 61-62: Set the anti-parallel rotation angle parameters; implemented by `sin_ang=0; cos_ang=-1`.
- Lines 67-68: Use the axis norm as the sine of the rotation angle; implemented by `sin_ang=axis_norm`.
- Lines 72-73: Normalise the rotation axis; implemented by `rot_axis=rot_axis/axis_norm`.
- Lines 75-76: Build the cross-product matrix of the rotation axis; implemented by `skew_mat=[0 -rot_axis(3) rot_axis(2)`.
- Lines 80-81: Build the rotation matrix using Rodrigues formula; implemented by `rot_mat=eye(3)+sin_ang*skew_mat+(1-cos_ang)*(skew_mat*skew_mat)`.

### Control flow inferred from the code

- Line 47: conditional branch on `axis_norm<align_tol`.
- Line 48: conditional branch on `cos_ang>0`.

### Key state/data transformations

- Lines 35: computes `v_from` using `v_from=v_from(:)/norm(v_from,2)`.
- Lines 36: computes `v_to` using `v_to=v_to(:)/norm(v_to,2)`.
- Lines 39: computes `align_tol` using `align_tol=1e-12`.
- Lines 42: computes `rot_axis` using `rot_axis=cross(v_from,v_to)`.
- Lines 43: computes `axis_norm` using `axis_norm=norm(rot_axis,2)`.
- Lines 44: computes `cos_ang` using `cos_ang=dot(v_from,v_to)`.
- Lines 51: computes `rot_mat` using `rot_mat=eye(3)`.
- Lines 57: computes `null_basis` using `null_basis=null(v_from.')`.
- Lines 62: computes `sin_ang` using `sin_ang=0; cos_ang=-1`.
- Lines 76: computes `skew_mat` using `skew_mat=[0 -rot_axis(3) rot_axis(2)`.

### Local helper functions

- Line 86: `grumble()` — `function grumble(v_from,v_to)`.
  - Representative operation: `if (~isnumeric(v_from))||(~isreal(v_from))||(~isvector(v_from))|| (numel(v_from)~=3)||any(~isfinite(v_from(:)))`.
  - Representative operation: `(numel(v_from)~=3)||any(~isfinite(v_from(:)))`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `v_from()`, `v_to()`, `cross()`, `dot()`, `null()`, `null_basis()`, `rot_axis()`, `isvector()`, `any()`, `eps()`.
