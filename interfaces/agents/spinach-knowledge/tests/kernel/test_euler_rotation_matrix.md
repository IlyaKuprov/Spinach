# tests/kernel/test_euler_rotation_matrix.m

- Signature: `result=test_euler_rotation_matrix()`

## Purpose

Tests active ZYZ Euler rotation matrices. Syntax: result=test_euler_rotation_matrix()

## Physical / mathematical content

## Numerical / algorithmic content

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
