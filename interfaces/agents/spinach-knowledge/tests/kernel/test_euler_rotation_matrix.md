# tests/kernel/test_euler_rotation_matrix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_euler_rotation_matrix.m`
- Signature: `result=test_euler_rotation_matrix()`
- Total lines: 39

## Purpose

Tests active ZYZ Euler rotation matrices. Syntax: result=test_euler_rotation_matrix()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Euler active rotation matrix\n')`.
- Lines 19-22: State the physical target of the test; implemented by `result=new_test_result('kernel/euler_rotation_matrix', 'Euler active rotation matrix', 'euler2dcm must implement the active ZYZ convention.')`.
- Lines 24-25: Build a simple ninety-degree Z rotation; implemented by `R=euler2dcm(pi/2,0,0)`.
- Lines 28-30: Check the active rotation and orthogonality; implemented by `result=test_close(result,'active Z rotation',R,R_ref,1e-15,1e-15, 'a positive active Z rotation maps the x axis into y')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/euler_rotation_matrix', 'Euler active rotation matrix', 'euler2dcm must implement the active ZYZ convention.')`.
- Lines 25: computes `R` using `R=euler2dcm(pi/2,0,0)`.
- Lines 26: computes `R_ref` using `R_ref=[0 -1 0;1 0 0;0 0 1]`.
- Lines 36: computes `'the documented action is v` using `'the documented action is v=R*v for column vectors')`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the active convention used by Spinach: alpha=pi/2,
- beta=0, gamma=0 is a counter-clockwise rotation around Z, taking x into y.

## Implementation structure

- Tests active ZYZ Euler rotation matrices. Syntax:
- result=test_euler_rotation_matrix()
- result -regression test result with explanatory messages
- The test checks the active convention used by Spinach: alpha=pi/2,
- beta=0, gamma=0 is a counter-clockwise rotation around Z, taking x into y.
- Announce the test target
- State the physical target of the test
- Build a simple ninety-degree Z rotation
- Check the active rotation and orthogonality

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `euler2dcm()`, `test_close()`.
