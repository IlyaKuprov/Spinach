# tests/kernel/test_transform_tensor_suite.m

- Signature: `result=test_transform_tensor_suite()`

## Purpose

Tests tensor transform helpers.

## Physical / mathematical content

- Checks the Haeberlen anisotropy and asymmetry, Mehring axiality and rhombicity, Herzfeld-Berger span and skew, and zero-field splitting tensor conventions.
- Checks quadrupolar tensor construction and conversion, electric-field-gradient scaling, rotational averaging, and spin-half Hamiltonian decomposition.

## Numerical / algorithmic content

- Tests interaction tensor parametrisations and round-trips between Cartesian and irreducible spherical tensor components.
- Checks isotropic–antisymmetric–symmetric decomposition and reconstruction, spherical harmonic coefficients of an isotropic quadratic form, and extraction of traceless symmetric matrix parameters.
- Compares computed values with reference tensors and reconstructions using the tolerances specified in the tests.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the test target and creates a result for `kernel/transform_tensor_suite`.
- Checks principal values for `anas2mat`, `axrh2mat`, `spsk2mat`, and `zfs2mat`, including zero-field splitting tracelessness.
- Checks `mat2axrh` values and eigenvalue order, then tests `mat2ias`/`ias2mat` and `mat2sphten`/`sphten2mat` round-trips.
- Checks `qform2sph`, `stev2sph`, and `tsm2param` against reference values and orientation reconstruction.
- Checks `eeqq2nqi`, `castep2nqi`, and both two-site `weblab2nqi` tensors.
- Checks `axis_tsymm` averaging around the z axis and `ham2nqi` decomposition of a spin-half Zeeman Hamiltonian.
