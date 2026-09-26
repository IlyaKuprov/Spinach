# tests/kernel/test_transform_rotation_suite.m

- Signature: `result=test_transform_rotation_suite()`

## Purpose

Tests rotation transform helpers. Syntax: result=test_transform_rotation_suite()

## Physical / mathematical content

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
