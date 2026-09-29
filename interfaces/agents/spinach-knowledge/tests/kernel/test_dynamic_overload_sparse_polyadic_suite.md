# tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m

## Purpose

Regression test for dynamic dispatch of sparse, polyadic, and OPIUM operator overloads in Spinach. The suite verifies that object-oriented operations on `rcv`, `polyadic`, and `opium` objects match explicit dense references computed on small deterministic operands.

Source: [tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m)

## Behaviour

The function announces the test target with `fprintf`, initialises a test result via `new_test_result` under the identifier `kernel/dynamic_overload_sparse_polyadic_suite`, and then runs assertion groups:

- **RCV construction and conversion**: builds deterministic sparse matrices `A` and `B`, wraps them as `rcv` objects, and checks `sparse(ra)`, `full(ra)`, `size(ra)`, and `size(ra,1)` against the source matrices with tolerances `1e-15` (exact match, tolerance zero, for size checks).
- **RCV arithmetic**: checks `ra+rb`, `ra+B`, `ra-rb`, `3.*ra`, `ra*2`, `2*ra`, and `ra./2` against dense references with tolerances `1e-15`.
- **RCV direct method dispatch**: calls `plus`, `minus`, `times`, `mtimes`, and `rdivide` directly by name to confirm dispatch to the RCV overloads.
- **RCV matrix products and transposition**: checks `ra*C`, `L*ra`, `ra.'`, conjugate transpose of a complex RCV built from `A+1i*B`, horizontal and vertical concatenation, direct `ctranspose`/`horzcat`/`vertcat`, and `gather` as a CPU no-op.
- **RCV plotting**: temporarily sets the default figure visibility to `off` (restored via `onCleanup`), calls `spy(ra)`, and asserts that at least one new figure was created before closing all figures.
- **Guarded RCV GPU path**: if `local_have_gpu` returns true, uploads `ra` with `gpuArray`, verifies the `isGPU` flag, and checks a `gather` round-trip; otherwise appends a `SKIP` message.
- **Polyadic construction and structure**: builds `polyadic` objects `P` and `Q` from 2-by-2 factor cell arrays, computes dense references `p_ref` and `q_ref` as sums of `kron` products, runs `validate(P)`, and checks `full(P)`, `size(P)`, and the predicates `isempty`, `allfinite`, and `nnz` (expected `nnz(P)` is 12).
- **Polyadic arithmetic and Kronecker dispatch**: checks `P+Q`, `P-Q`, `2*P`, `P*3`, `P*vec` for a dense vector, direct `plus`/`minus`/`mtimes`, sparse prefix/suffix multiplication with `sparse(eye(4))`, and `kron` with numeric left and right operands.
- **Polyadic transpose and structural operations**: checks `P.'`, conjugate transpose of a complex polyadic, direct `ctranspose`, `prefix(L,P)`, `suffix(P,R)`, `simplify` of identity-wrapped objects using `speye(4)`, and `inflate(P)`.
- **Guarded polyadic GPU path**: if a GPU is available, checks `gather(full(gpuArray(P)))` against `p_ref`; otherwise appends a `SKIP` message.
- **OPIUM operations**: builds `opium(3,1)`, `opium(3,2)`, and `opium(3,-4)`, and checks `size`, `full` of the unit operator against `eye(3)`, `sparse` conversion including the stored coefficient, left/right scalar multiplication, dense left/right multiplication with a 3-by-3 matrix `M`, OPIUM-by-OPIUM multiplication of coefficients, `kron` of two OPIUM objects (dimensions and coefficients multiply), `kron` with a numeric operand, and the predicates `allfinite`, `isnumeric`, `ismatrix`, and `nnz` (expected `nnz(scaled_op)` is 1).

The helper `local_have_gpu` returns false when `gpuDevice` is unavailable (checked with `exist('gpuDevice','file')~=2`) and otherwise attempts `gpuDevice()` inside a `try`/`catch`, returning true only on success.

## Inputs and outputs

- **Inputs**: none. The function takes no arguments.
- **Outputs**: `result` — regression test result structure with explanatory messages, produced by `new_test_result` and accumulated through `test_close` and `test_true` calls, including `PASS` and `SKIP` messages appended directly.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_sparse_polyadic_suite.m)
