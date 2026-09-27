# examples/fundamentals/tensor_structures/polyadic_test_1.m

- Signature: `polyadic_test_1()`

## Purpose

Checks that a nested `polyadic` matrix represents the same operator as its explicit dense Kronecker-product reference after prefix and suffix multiplication. It then exercises the object's arithmetic, matrix-vector operations, transpose/conjugate-transpose, Kronecker products, and Spinach `step` operation against dense results.

## Physical / mathematical content

- The test normalises complex matrix factors and vectors, builds a polyadic expression with nested and flat Kronecker terms, and applies random prefix and suffix matrices. The resulting dense reference is `M=2*p1*p2*kron(a,kron(b,c))*s1*s2`.
- This is an operator-algebra test rather than a specified spin-system simulation. Its final propagation check bootstraps Spinach and compares one `step` applied to the polyadic object and the reference matrix.

## Numerical / algorithmic content

- Dense equivalence checks use a one-norm difference threshold of `1e-14`; the `step` comparison uses `1e-12`.
- The test checks both left and right vector multiplication, additions and scalar arithmetic, conjugate-transpose identities, size preservation, and Kronecker-product identities.

## Implementation structure

- Generate and normalise complex factors, construct the nested polyadic, and apply its prefixes and suffixes.
- Compare `inflate` and `full` with the explicit dense matrix, then run the matrix-vector, addition, transpose, size, Kronecker, and `step` checks.
