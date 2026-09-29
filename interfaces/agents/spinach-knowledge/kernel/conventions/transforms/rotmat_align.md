# kernel/conventions/transforms/rotmat_align.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/rotmat_align.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rotmat_align.m)

- Signature: `rot_mat=rotmat_align(v_from,v_to)`

## Purpose

Return a `3x3` rotation matrix mapping the normalised `v_from` direction onto the normalised `v_to` direction. Inputs are directions; no angle or unit conversion is involved.

## Inputs and checks

- `v_from` and `v_to` must each be numeric, real, finite vectors with three elements. Row and column vectors are accepted.
- A vector whose 2-norm is below `eps('double')` is rejected. Both vectors are normalised before the alignment is calculated.

## Rotation convention and algorithm

The implementation forms `rot_axis=cross(v_from,v_to)`, `axis_norm=norm(rot_axis,2)`, and `cos_ang=dot(v_from,v_to)` from the normalised vectors. Its collinearity threshold is `1e-12`.

- If `axis_norm<1e-12` and `cos_ang>0`, the result is `eye(3)`.
- If `axis_norm<1e-12` and `cos_ang<=0`, it takes the first basis vector from `null(v_from')` as the anti-parallel rotation axis, with `sin_ang=0` and `cos_ang=-1`.
- Otherwise, the normalised cross product is the axis, `sin_ang=axis_norm`, and `cos_ang` is the dot product.

For the selected unit axis `a=(a1,a2,a3)`, the source builds

```text
S = [ 0  -a3   a2;
      a3    0  -a1;
     -a2   a1    0 ]
rot_mat = eye(3) + sin_ang*S + (1-cos_ang)*(S*S)
```

This is Rodrigues rotation. The source documents a minimum-angle alignment without an extra twist; in the anti-parallel case the axis is non-unique and the implementation selects one null-space basis vector. Output shape is `3x3`.
