# tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m

- Signature: `result=test_dynamic_overload_sparse_polyadic_suite()`

## Purpose

Tests dynamic dispatch of sparse, polyadic, and OPIUM overloads. Syntax: result=test_dynamic_overload_sparse_polyadic_suite()

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Outputs

- result -regression test result with explanatory messages
- The test exercises object-operation dispatch for rcv, polyadic,
- and opium objects using small deterministic dense references.

## Implementation structure

- Tests dynamic dispatch of sparse, polyadic, and OPIUM overloads. Syntax:
- result=test_dynamic_overload_sparse_polyadic_suite()
- result -regression test result with explanatory messages
- The test exercises object-operation dispatch for rcv, polyadic,
- and opium objects using small deterministic dense references.
- Announce the test target
- State the overload target of the test
- Build deterministic sparse operands
- Exercise RCV construction, conversion, and size dispatch
- Exercise RCV arithmetic dispatch
- Exercise direct RCV method-name dispatch for coverage tracking
- Exercise RCV matrix product and transposition dispatch
