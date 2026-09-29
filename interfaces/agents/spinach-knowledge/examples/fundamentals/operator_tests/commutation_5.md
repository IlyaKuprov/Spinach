# examples/fundamentals/operator_tests/commutation_5.m

- Signature: `commutation_5()`
- Source: [examples/fundamentals/operator_tests/commutation_5.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_5.m)

## Purpose

Checks finite-dimensional bosonic product and commutation relations, including the number-operator action on bosonic monomials.

## Operators and assumptions

The no-argument function chooses `nlevels=5+randi(20)`, so the cutoff is a random integer from 6 through 25 on each call. It obtains the Weyl operators from `weyl(nlevels)`, including annihilation `A.a`, creation `A.c`, number `A.n`, and `A.u`. It checks `A.c*A.a=A.n`; the number-operator commutators `[A.n,A.c]=A.c` and `[A.n,A.a]=-A.a`; and a boundary-adjusted commutator formed as `comm(A.a,A.c)` with its final diagonal entry set to one, compared with `A.u`. It then obtains `B=boson_mono(nlevels)`; for each monomial, `lin2kq(nlevels,n,1)` supplies `k,q` and the test checks `[A.n,B{n}]=(k-q)*B{n}`.

## Calling and numerical checks

Call `commutation_5()` in MATLAB with the cited Spinach functions available. Relative Frobenius-norm errors are compared with `acc=10*nlevels*eps('double')`. The source reports separate success messages for the Weyl product, the number-operator and boundary commutators, and the monomial commutators; an error is raised for a group whose error exceeds `acc`. These are conditional source-defined outcomes, not a report that the function was run.

The checks concern the finite matrices produced at one randomly selected cutoff per call. In particular, the explicit final-diagonal correction is part of the boundary test; the source does not test an untruncated oscillator or a range of user-specified cutoffs in one call.
