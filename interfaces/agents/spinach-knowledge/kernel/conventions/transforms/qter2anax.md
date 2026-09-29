# kernel/conventions/transforms/qter2anax.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/qter2anax.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=qter2anax.m)

- Signature: `[rot_axis,rot_angle]=qter2anax(q)`

## Behaviour

The function normalises `[q.u q.i q.j q.k]` by its Euclidean norm. Let `v=[q.i q.j q.k]` after normalisation and `n=norm(v,2)`. If `n==0`, it returns `rot_angle=0` and `rot_axis=[0 0 1]`. Otherwise it returns `rot_angle=2*atan2(n,q.u)` in radians and `rot_axis=v/n`.

## Input and outputs

The required structure fields are `u`, `i`, `j`, and `k`. The explicit checks reject missing fields and non-real concatenated components; this function does not explicitly check that the components are numeric scalars or that the quaternion norm is nonzero before normalisation.

- `rot_axis` — 1x3 row vector.
- `rot_angle` — scalar angle in radians.
