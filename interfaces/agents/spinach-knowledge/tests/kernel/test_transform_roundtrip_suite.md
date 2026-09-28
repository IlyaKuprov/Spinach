# tests/kernel/test_transform_roundtrip_suite.m

- Signature: `result=test_transform_roundtrip_suite()`

## Purpose

Tests deterministic coordinate and tensor transforms using exact geometrical identities, algebraic inverses, and known tensor decompositions.

## Physical / mathematical content

- Direction-cosine matrices are orthogonal proper rotations with determinant +1.
- Quaternion and angle-axis representations describe the same rotation. Euler angles may be ill-conditioned, but converting them back to a direction-cosine matrix reconstructs the active ZYZ rotation.
- With zero Euler angles, axiality/rhombicity conversion places the Mehring-order principal values on the diagonal. The isotropic part is their mean; axiality is `2*zz-(xx+yy)` and rhombicity is `yy-xx`.
- The nine irreducible spherical tensor components span all `3x3` Cartesian tensors.
- Fractional coordinates in an orthorhombic cell scale by the cell edges. Spherical coordinates use the ISO radius/inclination/azimuth convention: inclination is measured from positive `z`, and azimuth in the `xy` plane from positive `x`.

## Numerical / algorithmic content

- Checks direction-cosine-matrix orthogonality and determinant, quaternion/angle-axis rotation round-tripping, and Euler-to-matrix reconstruction with absolute and relative tolerances of `1e-14`.
- Checks `axrh2mat` principal values and `mat2axrh` isotropic part, axiality, and rhombicity with tolerances of `1e-14`.
- Converts a general `3x3` Cartesian tensor through `mat2sphten` and `sphten2mat`, checking reconstruction with tolerances of `1e-13`.
- Checks `frac2cart` coordinates and primitive vectors for an orthorhombic cell with edges `2`, `3`, and `4` and angles `90`, `90`, and `90` degrees, using tolerances of `1e-14`.
- Checks `xyz2sph` radii, inclinations, and azimuths for the Cartesian basis vectors with tolerances of `1e-14`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the coordinate and tensor transform test target and initializes a regression test result.
- Runs rotation, tensor-decomposition, fractional-coordinate, and spherical-coordinate checks through `test_close`.
