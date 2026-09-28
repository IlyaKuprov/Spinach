# examples/fundamentals/convention_tests/stevens_test.m

- Signature: `stevens_test()`

## Purpose

Tests Spinach's Stevens-operator matrices against explicit expressions from the literature, and checks the rank-2 operators against irreducible spherical tensor (IST) expressions.

## Method and checks

With spin quantum number 12 (multiplicity 25), the test compares every component q=−k,…,k for ranks k=6, 4, and 2, using explicit expressions built from the spin operators. The 1-norm of the rank-wise residual vector must be at most 10⁻⁴, 10⁻⁷, and 10⁻¹², respectively. It then compares the five rank-2 Stevens matrices with the corresponding IST combinations; the residual-vector 1-norm must be at most 10⁻¹⁰.
