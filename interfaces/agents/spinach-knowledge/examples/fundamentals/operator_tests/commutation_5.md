# examples/fundamentals/operator_tests/commutation_5.m

- Signature: `commutation_5()`

## Purpose

Tests finite-dimensional bosonic product and commutation relations.

## Physical / mathematical content

A random truncated-boson dimension is chosen as `5+randi(20)`. From the Weyl operators, the test checks the number-operator product `c*a=n`, the number-operator commutators `[n,c]=c` and `[n,a]=-a`, and the finite-boundary commutator after setting its final diagonal entry to one. It then checks `[n,B]=(k-q)B` for every bosonic monomial returned by `boson_mono`, with indices obtained from `lin2kq`.

## Numerical / algorithmic content

The accuracy threshold is `10*nlevels*eps('double')`. Product and commutator discrepancies are normalized by the Frobenius norm of the reference operator; exceeding the threshold raises an error.

## Implementation structure

The example performs three groups of checks: the Weyl product, number-operator and corrected boundary commutators, and number-operator commutators with all bosonic monomials. Each group reports success or stops with its corresponding failure.
