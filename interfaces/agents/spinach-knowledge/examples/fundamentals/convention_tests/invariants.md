# examples/fundamentals/convention_tests/invariants.m

- Signature: `invariants()`

## Purpose

Tests the invariant identity associated with Equation 3 in [doi:10.1002/chem.200902300](http://dx.doi.org/10.1002/chem.200902300), for randomly generated rhombic and axial eigenvalue sets.

## Checks

For 100 random symmetric matrices, the code forms axiality `Ax=2*D3-(D1+D2)` and rhombicity `Rh=D1-D2` from shuffled eigenvalues. It compares `(Ax^2+3*Rh^2)/6` with the equivalent eigenvalue expression; the absolute difference must be at most 10⁻⁶.

A second 100-case loop constructs zero-trace axial eigenvalues, shuffles them, and compares `(D1-D2)^2` with `D1^2+D2^2+D3^2-D1*D2-D1*D3-D2*D3`, again with an absolute tolerance of 10⁻⁶.
