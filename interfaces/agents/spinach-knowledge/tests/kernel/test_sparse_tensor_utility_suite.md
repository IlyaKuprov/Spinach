# tests/kernel/test_sparse_tensor_utility_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_sparse_tensor_utility_suite.m`
- Signature: `result=test_sparse_tensor_utility_suite()`
- Total lines: 69

## Purpose

Tests sparse, tensor, and numerical utility helpers. Syntax: result=test_sparse_tensor_utility_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `kronm()`, `kronm_new()`, `blinv()`, `blprod()`, `spden()`, `frob_chop()`, `svd_shrink()`.
