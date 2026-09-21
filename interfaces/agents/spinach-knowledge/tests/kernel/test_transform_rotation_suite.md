# tests/kernel/test_transform_rotation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_rotation_suite.m`
- Signature: `result=test_transform_rotation_suite()`
- Total lines: 138

## Purpose

Tests rotation transform helpers. Syntax: result=test_transform_rotation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks active ZYZ Euler rotations, Euler/DCM inversion, rotation
- composition, angle-axis normalisation, quaternion round-trips, Wigner
- identity rotation, and minimum-angle vector alignment.

## Implementation structure

- Tests rotation transform helpers. Syntax:
- result=test_transform_rotation_suite()
- result -regression test result with explanatory messages
- The test checks active ZYZ Euler rotations, Euler/DCM inversion, rotation
- composition, angle-axis normalisation, quaternion round-trips, Wigner
- identity rotation, and minimum-angle vector alignment.
- Announce the test target
- State the geometric target of the test
- Check the active ZYZ convention on a known quarter-turn
- Compare equivalent rotations through their Euler angle degeneracy
- Check angular tolerance acceptance
- Check angular tolerance rejection

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `euler2dcm()`, `test_close()`, `test_true()`, `euler_equiv()`, `dcm2euler()`, `euler_sup()`, `anax2dcm()`, `anax2qter()`, `qter2anax()`, `euler2qter()`, `qter2dcm()`, `qter2euler()`, `dcm2qter()`, `dcm2wigner()`, `rotmat_align()`.
