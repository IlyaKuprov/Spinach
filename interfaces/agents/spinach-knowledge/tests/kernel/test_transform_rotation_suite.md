# tests/kernel/test_transform_rotation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_rotation_suite.m`
- Signature: `result=test_transform_rotation_suite()`
- Total lines: 138

## Purpose

Tests rotation transform helpers. Syntax: result=test_transform_rotation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Rotation transform helpers\n')`.
- Lines 20-23: State the geometric target of the test; implemented by `result=new_test_result('kernel/transform_rotation_suite', 'Rotation transform helpers', 'rotation transforms must preserve active rotation, composition, and round-trip c…`.
- Lines 25-26: Check the active ZYZ convention on a known quarter-turn; implemented by `R=euler2dcm(pi/2,0,0)`.
- Lines 33-35: Compare equivalent rotations through their Euler angle degeneracy; implemented by `result=test_true(result,'euler_equiv zero-beta degeneracy',euler_equiv([0.3 0 0.4],[0.7 0 0],1e-14), 'for beta equal to zero, alpha and gamma only enter through their su…`.
- Lines 37-39: Check angular tolerance acceptance; implemented by `result=test_true(result,'euler_equiv tolerance pass',euler_equiv([0 0 0],[5e-7 0 0],1e-6), 'relative rotations inside the angular tolerance must pass')`.
- Lines 41-43: Check angular tolerance rejection; implemented by `result=test_true(result,'euler_equiv tolerance fail',~euler_equiv([0 0 0],[2e-6 0 0],1e-6), 'relative rotations outside the angular tolerance must fail')`.
- Lines 45-46: Recover a non-singular Euler rotation through its DCM; implemented by `angles=[0.37 0.91 -0.42]`.
- Lines 59-60: Compose two active rotations in the documented order; implemented by `rot_one=[0.2 0.4 -0.3]`.
- Lines 68-69: Check angle-axis normalisation and inverse rotation; implemented by `axis_vec=[1;2;3]`.
- Lines 87-88: Round-trip a non-zero angle through quaternion form; implemented by `q=anax2qter(axis_vec,angle)`.
- Lines 99-100: Check the quaternion conversions against the Euler convention; implemented by `qtr=euler2qter(0.4,1.1,-0.7)`.
- Lines 111-112: Check that identity rotation has identity second-rank Wigner matrix; implemented by `D=dcm2wigner(eye(3))`.
- Lines 118-119: Align two non-zero vectors and check the proper rotation properties; implemented by `v_from=[2;0;0]`.
- Lines 129-130: Check the anti-parallel branch explicitly; implemented by `R=rotmat_align([1;0;0],[-1;0;0])`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/transform_rotation_suite', 'Rotation transform helpers', 'rotation transforms must preserve active rotation, composition, and round-trip c…`.
- Lines 26: computes `R` using `R=euler2dcm(pi/2,0,0)`.
- Lines 27: computes `R_ref` using `R_ref=[0 -1 0;1 0 0;0 0 1]`.
- Lines 46: computes `angles` using `angles=[0.37 0.91 -0.42]`.
- Lines 48: computes `angles_obs` using `angles_obs=dcm2euler(R)`.
- Lines 52: computes `'the beta` using `'the beta=0 gimbal case must reconstruct exactly')`.
- Lines 55: computes `R_dirty` using `R_dirty=R+1e-7*[0.3 -0.7 0.2; 0.1 0.4 -0.6; -0.5 0.2 0.3]`.
- Lines 60: computes `rot_one` using `rot_one=[0.2 0.4 -0.3]`.
- Lines 61: computes `rot_two` using `rot_two=[-0.5 0.7 0.6]`.
- Lines 62: computes `rot_cmp` using `rot_cmp=euler_sup(rot_one,rot_two)`.
- Lines 63: computes `R_cmp` using `R_cmp=euler2dcm(rot_cmp)`.
- Lines 69: computes `axis_vec` using `axis_vec=[1;2;3]`.
- Lines 70: computes `angle` using `angle=0.73`.
- Lines 71: computes `R_one` using `R_one=anax2dcm(axis_vec,angle)`.
- Lines 72: computes `R_two` using `R_two=anax2dcm(5*axis_vec,angle)`.
- Lines 73: computes `R_inv` using `R_inv=anax2dcm(axis_vec,-angle)`.
- Lines 88: computes `q` using `q=anax2qter(axis_vec,angle)`.
- Lines 89: computes `[axis_obs,angle_obs]` using `[axis_obs,angle_obs]=qter2anax(q)`.

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
