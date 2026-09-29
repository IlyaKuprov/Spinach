# tests/kernel/test_transform_roundtrip_suite.m

## Purpose

Regression test suite for Spinach coordinate and tensor transformation functions. The suite verifies that transformation helpers preserve rotations, coordinates, and tensor decompositions, using exact geometrical identities, algebraic inverses, and known tensor decompositions.

## Behaviour

The function announces the test target with `fprintf`, then initialises a regression test result object via `new_test_result` with the identifier `kernel/transform_roundtrip_suite`, the description "Coordinate and tensor transform functions", and the requirement that "transform helpers must preserve rotations, coordinates, and tensor decompositions". Each check is registered through `test_close`, which compares computed values against references at tight tolerances (mostly `1e-14`, spherical tensor round-trip at `1e-13`) and appends explanatory messages to the result.

Checks performed:

- **`anax2dcm` orthogonality and determinant**: the direction-cosine matrix from axis `[0 0 2]` and angle `pi/3` must satisfy `R'*R = eye(3)` and `det(R) = 1` (proper rotations have determinant +1).
- **`anax2qter`/`qter2anax` round-trip**: converting axis `[1 2 3]`, angle `0.37*pi` to a quaternion and back must reproduce the same rotation matrix as direct `anax2dcm` on the original axis and angle.
- **`dcm2euler`/`euler2dcm` round-trip**: for Euler angles `[0.21*pi 0.37*pi 0.43*pi]`, recovering angles from the DCM and reconstructing must return the original active ZYZ rotation matrix; the comment notes Euler conversion is ill-conditioned in angles but DCM reconstruction is unique.
- **`axrh2mat` principal values**: with isotropic value `4`, axiality `6`, rhombicity `2`, and zero Euler angles, the matrix must equal `diag([iso-(ax+3*rh)/6, iso-(ax-3*rh)/6, iso+ax/3])` (Mehring-order eigenvalues on the diagonal).
- **`mat2axrh` decomposition**: from the reference diagonal matrix, the recovered isotropic part must equal the mean principal value, axiality must equal `2*eigvals(3)-(eigvals(1)+eigvals(2))`, and rhombicity must equal `eigvals(2)-eigvals(1)`, in Mehring order.
- **`mat2sphten`/`sphten2mat` round-trip**: for the Cartesian tensor `[1 2 3;4 5 6;7 8 10]`, converting to irreducible spherical tensor ranks 0, 1, 2 and back must reproduce the original matrix; the nine spherical components span all 3x3 Cartesian tensors.
- **`frac2cart` orthorhombic cell**: with cell edges `2, 3, 4` and angles `90, 90, 90`, fractional coordinates `ABC=[0 0 0; 1/2 1/3 1/4; 1 1 1]` must map to `ABC*diag([2 3 4])`, and the primitive vectors must be `diag([2 3 4])`.
- **`xyz2sph` ISO convention**: for the Cartesian unit basis vectors, radius must be `[1 1 1]`, inclination `[pi/2 pi/2 0]` (measured down from the positive z axis), and azimuth `[0 pi/2 0]` (measured in the xy plane from the positive x axis).

## Inputs and outputs

- **Outputs**: `result` — regression test result object with explanatory messages, produced by `new_test_result` and accumulated `test_close` comparisons.
- **Inputs**: none.

## References

- Source: [tests/kernel/test_transform_roundtrip_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_roundtrip_suite.m)
