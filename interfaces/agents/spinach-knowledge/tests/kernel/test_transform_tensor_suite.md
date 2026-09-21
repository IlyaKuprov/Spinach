# tests/kernel/test_transform_tensor_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_tensor_suite.m`
- Signature: `result=test_transform_tensor_suite()`
- Total lines: 149

## Purpose

Tests tensor transform helpers. Syntax: result=test_transform_tensor_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks interaction tensor parametrisations, spherical tensor
- round-trips, quadrupolar conversions, axial symmetrisation, and simple
- Hamiltonian decomposition.

## Implementation structure

- Tests tensor transform helpers. Syntax:
- result=test_transform_tensor_suite()
- result -regression test result with explanatory messages
- The test checks interaction tensor parametrisations, spherical tensor
- round-trips, quadrupolar conversions, axial symmetrisation, and simple
- Hamiltonian decomposition.
- Announce the test target
- State the tensor target of the test
- Check Haeberlen anisotropy and asymmetry at zero Euler angles
- Check axiality and rhombicity matrix construction
- Check span and skew construction in the Herzfeld-Berger convention
- Check zero-field splitting tensor construction

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `anas2mat()`, `test_close()`, `axrh2mat()`, `mat2axrh()`, `spsk2mat()`, `zfs2mat()`, `mat2ias()`, `ias2mat()`, `mat2sphten()`, `sphten2mat()`, `qform2sph()`, `stev2sph()`, `tsm2param()`, `euler2dcm()`, `eeqq2nqi()`.
