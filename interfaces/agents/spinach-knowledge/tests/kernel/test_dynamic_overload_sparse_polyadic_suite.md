# tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m`
- Signature: `result=test_dynamic_overload_sparse_polyadic_suite()`
- Total lines: 253

## Purpose

Tests dynamic dispatch of sparse, polyadic, and OPIUM overloads. Syntax: result=test_dynamic_overload_sparse_polyadic_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file also defines local helper function(s): `local_have_gpu()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `rcv()`, `test_close()`, `double()`, `plus()`, `minus()`, `times()`, `mtimes()`, `rdivide()`, `ctranspose()`, `horzcat()`, `vertcat()`, `gather()`, `get()`, `set()`, `onCleanup()`.
