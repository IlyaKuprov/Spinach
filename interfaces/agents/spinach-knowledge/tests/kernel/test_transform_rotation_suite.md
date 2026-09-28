# tests/kernel/test_transform_rotation_suite.m

- Signature: `result=test_transform_rotation_suite()`

## Purpose

Tests rotation transform helpers.

## Physical / mathematical content

The test checks active ZYZ Euler rotations, Euler/DCM inversion, rotation composition, angle-axis normalisation, quaternion round-trips, Wigner identity rotation, and minimum-angle vector alignment.

## Numerical / algorithmic content

Checks rotations against reference matrices and round-trip reconstructions using the tolerances specified in the source tests.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Checks the active ZYZ convention with a quarter-turn and verifies equivalent Euler rotations and angular tolerance behavior.
- Reconstructs rotations from DCMs, including gimbal cases and a subtly corrupted input, and checks rotation composition.
- Tests angle-axis normalisation, inverse rotations, orthogonality, determinant, and agreement with active Euler rotations.
- Tests quaternion round-trips and agreement with Euler and DCM representations.
- Checks the second-rank Wigner matrix for identity rotation and its unitarity.
- Tests minimum-angle vector alignment, including the anti-parallel case.
