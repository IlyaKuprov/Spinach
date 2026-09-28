# tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m

- Signature: `result=test_dynamic_overload_sparse_polyadic_suite()`

## Purpose

Regression test for dynamic dispatch of `rcv`, `polyadic`, and `opium` overloads against small deterministic dense references. Returns a test result with explanatory messages.

## Numerical / algorithmic content

- Constructs `rcv` objects from deterministic sparse matrices and checks conversion, dimensions, arithmetic, direct method-name dispatch, matrix products, transpose and conjugate transpose, concatenation, and CPU `gather` against matrix references. It also checks that `spy` creates a figure with plotting set offscreen.
- Constructs deterministic `polyadic` objects and compares their opened values with explicit sums of Kronecker products. Checks validation, dimensions, structural predicates, arithmetic, vector multiplication, direct dispatch, sparse prefix and suffix multiplication, Kronecker products, transpose and conjugate transpose, `simplify`, and `inflate`.
- Checks `opium` dimensions, full and sparse conversion, scalar and matrix multiplication, Kronecker products, and structural predicates against identity-matrix references.
- Uses absolute and relative tolerances of `1e-15` for value comparisons and zero tolerances for dimension comparisons. GPU upload and `gather` checks for `rcv` and `polyadic` run only when a usable GPU is available; otherwise, the result records skip messages.