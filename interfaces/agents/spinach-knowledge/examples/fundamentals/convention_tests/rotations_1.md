# examples/fundamentals/convention_tests/rotations_1.m

- MATLAB implementation: [examples/fundamentals/convention_tests/rotations_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/rotations_1.m)

- Signature: `rotations_1()`

## Purpose

A seven-check internal consistency suite for Spinach's rotation and Cartesian/spherical-tensor conversion routines. It distinguishes DCM/Euler/Wigner conventions from tensor round trips and quaternion/angle-axis conversions, which helps localise which representation boundary a mismatch involves.

## Setup and checks

Run `rotations_1()`. One random symmetric traceless 3×3 matrix `A` and one Euler triple (drawn componentwise as `rand(1,3).*[2*pi pi 2*pi]`) are reused across the checks:

1. Convert `euler2dcm(eulers)*A*euler2dcm(eulers)'` to rank-2 spherical components and compare with rank-2 components of `A` rotated by `wigner(2,...)`; the 2-norm residual must be below `1e-10`.
2. Convert Euler angles to a DCM and back with `dcm2euler`; the direct angle-vector difference must have 2-norm below `1e-3`.
3. Compare `dcm2wigner(euler2dcm(eulers))` with `wigner(2,...)`; 2-norm residual below `1e-10`.
4. Round-trip `A` through `mat2sphten` and `sphten2mat`; matrix residual below `1e-10` in the 2-norm.
5. Convert the Euler DCM to a quaternion with `dcm2qter` and back with `qter2dcm`; DCM residual below `1e-10`.
6. Normalise a random quaternion, convert it to angle-axis with `qter2anax` and back with `anax2qter`; the four quaternion-component residual has 2-norm below `1e-10`.
7. For a normalised random quaternion, compare its direct `qter2dcm` DCM with the DCM made from its `qter2anax` angle-axis pair using `anax2dcm`; residual below `1e-10`.

## Observable result and scope

Each successful check prints `Test 1 passed.` through `Test 7 passed.`; the first failed check raises its numbered inconsistency error. The function produces no plot and returns no declared result. It draws a new matrix, Euler triple, and quaternions per run without setting a seed; the result is sampled internal-consistency coverage, not exhaustive validation of every rotation or singular convention.
