# tests/kernel/test_transform_tensor_suite.m

- Signature: `result=test_transform_tensor_suite()`

## Purpose

Tests tensor transform helpers. Syntax: result=test_transform_tensor_suite()

## Physical / mathematical content

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
