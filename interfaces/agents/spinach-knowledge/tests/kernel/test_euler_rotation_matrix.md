# tests/kernel/test_euler_rotation_matrix.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_euler_rotation_matrix.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_euler_rotation_matrix.m)

## Purpose

Regression test for the active ZYZ Euler rotation matrices produced by `euler2dcm`. The test verifies that Spinach's active convention is implemented correctly: `alpha=pi/2`, `beta=0`, `gamma=0` is a counter-clockwise rotation around Z, taking x into y.

## Behaviour

- Announces the test target with `fprintf('TESTING: Euler active rotation matrix\n')`.
- Initialises a regression test result via `new_test_result` for the target `kernel/euler_rotation_matrix`, with the description "Euler active rotation matrix" and the specification "euler2dcm must implement the active ZYZ convention."
- Builds a ninety-degree Z rotation with `R=euler2dcm(pi/2,0,0)` and the reference matrix `R_ref=[0 -1 0;1 0 0;0 0 1]`.
- Runs four checks with `test_close`, each using tolerances `1e-15` (absolute and relative):
  - **active Z rotation**: `R` against `R_ref`, with the message "a positive active Z rotation maps the x axis into y".
  - **orthogonality**: `R'*R` against `eye(3)`, with the message "proper rotations preserve vector lengths".
  - **proper determinant**: `det(R)` against `1`, with the message "direction cosine matrices must have determinant +1".
  - **vector action**: `R*[1;0;0]` against `[0;1;0]`, with the message "the documented action is v=R*v for column vectors".

## Inputs and outputs

- **Inputs:** none. The function is called as `result=test_euler_rotation_matrix()`.
- **Outputs:**
  - `result` — regression test result with explanatory messages.

## References

- `euler2dcm` — builds the direction cosine matrix under test.
- `new_test_result` — initialises the regression test result structure.
- `test_close` — performs the numerical comparisons and records messages.
