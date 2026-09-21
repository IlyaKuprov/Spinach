# tests/kernel/test_euler_rotation_matrix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_euler_rotation_matrix.m`
- Signature: `result=test_euler_rotation_matrix()`
- Total lines: 39

## Purpose

Tests active ZYZ Euler rotation matrices. Syntax: result=test_euler_rotation_matrix()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `euler2dcm()`, `test_close()`.
