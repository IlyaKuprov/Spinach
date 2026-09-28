# kernel/conventions/transforms/euler_sup.m

- Signature: `rot_cmp=euler_sup(rot_one,rot_two)`

## Purpose

Composes two ZYZ active Euler rotations. The rotations act in the supplied order: `v_rot = R_two * R_one * v`, so `R_comp = R_two * R_one`.

## Parameters / inputs

- `rot_one`, `rot_two`: real, finite three-element Euler-angle vectors `[alpha beta gamma]` in radians, using the ZYZ active convention.

## Output

- `rot_cmp`: row vector `[alpha beta gamma]` for the composite rotation, in radians. The angles are wrapped into `(-pi,pi]`.

## Method

The function builds both direction-cosine matrices and multiplies them as `dcm_two*dcm_one`, then converts the composite matrix back to Euler angles. Identity and the singular `beta=pi` branch are handled explicitly.

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=euler_sup.m)
