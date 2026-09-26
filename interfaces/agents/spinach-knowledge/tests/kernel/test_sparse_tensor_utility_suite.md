# tests/kernel/test_sparse_tensor_utility_suite.m

- Signature: `result=test_sparse_tensor_utility_suite()`

## Purpose

Tests sparse, tensor, and numerical utility helpers. Syntax: result=test_sparse_tensor_utility_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks Kronecker-product application, Blicharski invariants,
- spectral densities, SVD truncation helpers, sparse density, and related
- numerical utilities against direct definitions.

## Implementation structure

- Tests sparse, tensor, and numerical utility helpers. Syntax:
- result=test_sparse_tensor_utility_suite()
- result -regression test result with explanatory messages
- The test checks Kronecker-product application, Blicharski invariants,
- spectral densities, SVD truncation helpers, sparse density, and related
- numerical utilities against direct definitions.
- Announce the test target
- State the utility target of the test
- Kronecker-product application must match explicit kron products
- Blicharski invariants ignore isotropic trace and follow their definitions
- Lorentzian spectral density follows its closed definition
- Frobenius SVD truncation keeps the smallest rank whose dropped tail is below tolerance
