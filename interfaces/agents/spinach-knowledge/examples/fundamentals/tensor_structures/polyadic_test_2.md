# examples/fundamentals/tensor_structures/polyadic_test_2.m

- Signature: `polyadic_test_2()`

## Purpose

Validates the `polyadic` matrix object's constructor, dense conversion, and overloaded operations against explicit Kronecker-product matrix references. The test covers composition, arithmetic, sparse and dense operands, empty matrices, nested simplification, and GPU conversion when a device is available.

## Physical / mathematical content

- The main reference is the sum `kron(a,b)+kron(c,d)`, represented as a two-term polyadic object. Further references are formed by ordinary dense matrix operations on that sum.
- This is a matrix-algebra test; it does not define a physical spin system or dynamics.

## Numerical / algorithmic content

- Most dense-reference identities are checked at `1e-12`; the optional GPU round-trip uses `1e-10`. GPU coverage is skipped when `gpuDeviceCount` is zero.
- Checks include constructor/`full`/`inflate` consistency, `validate`, prefix and suffix multiplication, size and emptiness, addition and subtraction, scalar and matrix multiplication, Kronecker products, transpose operations, finiteness, nonzero counts, and simplification of nested expressions.

## Implementation structure

- Construct complex dense and sparse factors, form a two-term polyadic object, and compare it with the dense sum-of-Kronecker-products reference.
- Apply each operation to both polyadic and dense forms, asserting their results agree.
- Exercise the zero-dimension and nested-simplification cases, then test GPU upload only if hardware is available.
