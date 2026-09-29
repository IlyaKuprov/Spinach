# examples/fundamentals/tensor_structures/polyadic_test_1.m

- MATLAB implementation: [examples/fundamentals/tensor_structures/polyadic_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/tensor_structures/polyadic_test_1.m)

- Signature: `polyadic_test_1()`
- Source: [`examples/fundamentals/tensor_structures/polyadic_test_1.m`](../../../../../../examples/fundamentals/tensor_structures/polyadic_test_1.m)

## Purpose

Checks that a nested `polyadic` representation agrees with a dense Kronecker-product reference, then tests its matrix/vector operations, arithmetic, conjugate transpose, size behaviour, Kronecker products, and Spinach `step` action.

## Tensor construction and scope

The test uses complex factors `a` (7 by 7), `b` (1 by 1), and sparse `c` (9 by 9, density 0.1), plus a complex vector of length 63. The factors and vector are normalised, and the polyadic is nested from `a`, `b`, and `c`; random normalised prefixes and suffixes are then applied. Its dense reference is `M = 2*p1*p2*kron(a,kron(b,c))*s1*s2`. Tests compare `inflate(P)` and `full(P)` with `M`, and compare left and right matrix-vector products.

This is a tensor/operator representation check rather than a defined physical spin model: the source does not specify spins, interactions, Hamiltonian, or a basis. It bootstraps `spin_system=bootstrap('hush')` only for the final `step` comparison; no physical model settings or unit convention are supplied.

## Operations and numerical settings

The source checks addition, scaling, conjugate-transpose behaviour, equality of the reported and dense sizes, and Kronecker-product identities. A random complex 2 by 2 matrix is normalised for the Kronecker checks. These algebraic comparisons use one-norm tolerances of 1e-14. For the final step check it sets `dt=1/norm(M,1)` and requires the vector difference between `step(spin_system,P,v,dt)` and `step(spin_system,M,v,dt)` to be below 1e-12. `dt` is a numerical value derived from the matrix norm, not a documented physical time unit.

## Use and limits

Run `polyadic_test_1` in the Spinach MATLAB environment. It reports success or raises an error independently for each assertion. The thresholds and expected outcomes in the source are test criteria; they are not evidence that the test has been run or passed.
