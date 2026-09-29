# kernel/conventions/transforms/euler_sup.m

**MATLAB source:** [kernel/conventions/transforms/euler_sup.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/euler_sup.m)
**Spinach Wiki:** [euler_sup.m](https://spindynamics.org/wiki/index.php?title=euler_sup.m)

## Purpose and convention

Composes two ZYZ active Euler rotations. The supplied order is the order of action on a vector: `v_rot = R_two * R_one * v`, so the composite rotation matrix is `R_comp = R_two * R_one`. Angles and returned values are in radians.

## Inputs and output

Signature: `rot_cmp=euler_sup(rot_one,rot_two)`

- `rot_one` and `rot_two` are numeric, real, finite three-element vectors `[alpha beta gamma]` in the ZYZ active convention. Row and column vectors are accepted.
- `rot_cmp` is a row vector `[alpha beta gamma]` for the composite rotation.

## Matrix composition and angle recovery

Each input is reshaped to a row and wrapped with `wrapToPi`. The source builds `dcm_one=euler2dcm(rot_one)` and `dcm_two=euler2dcm(rot_two)`, then forms `dcm_comp=dcm_two*dcm_one`. In the general case, [dcm2euler.m](dcm2euler.md) recovers the angles. The output is then wrapped into `(-pi,pi]`.

The implementation has two explicit branches, with threshold `1e-12`:

- If every diagonal entry of `dcm_comp` differs from 1 by less than `1e-12`, it returns the identity representative `[pi/2 0 -pi/2]`.
- If `abs(dcm_comp(3,3)+1) < 1e-12` and the first and third components of the two wrapped input triples differ by less than `1e-12`, it uses `alpha_m_gamma=atan2(-dcm_comp(2,1),-dcm_comp(1,1))`, sets `gam=-alpha_m_gamma/2`, and returns `[-gam pi gam]` before the final wrapping.

Euler representations are not unique; these branches select explicit representatives while preserving the composed matrix.
