# tests/kernel/test_transform_rotation_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_rotation_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_rotation_suite.m)

## Purpose

Regression test suite for the rotation transform helpers in the kernel module. The suite verifies that rotation transforms preserve active rotation, composition, and round-trip conventions, covering active ZYZ Euler rotations, Euler/DCM inversion, rotation composition, angle-axis normalisation, quaternion round-trips, Wigner identity rotation, and minimum-angle vector alignment.

## Behaviour

- Announces the test target with `fprintf('TESTING: Rotation transform helpers\n')` and initialises a test result object via `new_test_result` for `kernel/transform_rotation_suite`.
- **euler2dcm:** checks the active ZYZ convention on a quarter turn, verifying that `euler2dcm(pi/2,0,0)` equals `[0 -1 0;1 0 0;0 0 1]` (a positive active Z rotation maps the x axis into y), and that the one-vector syntax `euler2dcm([pi/2 0 0])` matches the three-scalar syntax, both at tolerances `1e-15`.
- **euler_equiv:** verifies the zero-beta degeneracy (`[0.3 0 0.4]` vs `[0.7 0 0]` at `1e-14`), tolerance acceptance (`[0 0 0]` vs `[5e-7 0 0]` at `1e-6`), and tolerance rejection (`[0 0 0]` vs `[2e-6 0 0]` at `1e-6`).
- **dcm2euler:** reconstructs a non-singular Euler rotation `[0.37 0.91 -0.42]` through its DCM at `1e-14`; checks the beta=0 and beta=pi gimbal cases with `euler2dcm(0.4,0,1.1)` and `euler2dcm(0.4,pi,1.1)` reconstruct exactly at `1e-14`; and checks that a DCM corrupted by a `1e-7`-scaled perturbation matrix returns the angles of the nearest proper rotation at `1e-6`.
- **euler_sup:** composes two active rotations `rot_one=[0.2 0.4 -0.3]` and `rot_two=[-0.5 0.7 0.6]`, requiring the composite matrix to equal `euler2dcm(rot_two)*euler2dcm(rot_one)` for column vectors, at `1e-7` tolerance.
- **anax2dcm:** with axis `[1;2;3]` and angle `0.73`, verifies axis normalisation (scaling the axis by 5 does not change the rotation, `1e-15`), inverse rotation composing to identity (`1e-14`), proper orthogonality `R'*R=eye(3)` (`1e-14`), determinant +1 (`1e-14`), and active convention agreement with ZYZ Euler rotations about the Z axis (`anax2dcm([0 0 1],0.31)` vs `euler2dcm(0,0,0.31)`) and Y axis (`anax2dcm([0 1 0],0.27)` vs `euler2dcm(0,0.27,0)`), both at `1e-14`.
- **anax2qter / qter2anax:** round-trips the same axis-angle through quaternion form, checking unit quaternion norm (`1e-15`), recovered axis equal to the normalised input axis (`1e-15`), and recovered angle equal to the input angle for angles below pi without branch ambiguity (`1e-15`).
- **euler2qter / qter2dcm / qter2euler / dcm2qter:** for Euler angles `(0.4,1.1,-0.7)`, checks unit quaternion norm (`1e-15`), that `qter2dcm(qtr)` reproduces the active ZYZ rotation matrix (`1e-14`), that `qter2euler` recovers angles giving the same rotation (`1e-14`), and that the DCM round trip `qter2dcm(dcm2qter(...))` is exact to machine precision (`1e-14`).
- **dcm2wigner:** checks that the identity rotation yields the identity second-rank Wigner matrix `eye(5)` (`1e-15`) and that `D'*D=eye(5)` (unitarity, `1e-15`).
- **rotmat_align:** aligns `v_from=[2;0;0]` to `v_to=[0;3;0]`, verifying the matrix takes the normalised source vector into the normalised target vector (`1e-15`), orthogonality (`1e-15`), and determinant +1 so no reflection is introduced (`1e-15`). The anti-parallel branch is checked explicitly with `[1;0;0]` to `[-1;0;0]`, requiring a pi rotation around an orthogonal axis (`1e-15`) with determinant +1 (`1e-15`).

## Inputs and outputs

- **Syntax:** `result=test_transform_rotation_suite()`
- **Outputs:**
  - `result` — regression test result with explanatory messages.
- **Inputs:** none.

## References

- [Spinach source: tests/kernel/test_transform_rotation_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_transform_rotation_suite.m)
