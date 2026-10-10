# examples/fundamentals/convention_tests/invariants.m

- MATLAB implementation: [examples/fundamentals/convention_tests/invariants.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/invariants.m)

- Signature: `invariants()`

## Purpose

Checks two algebraic forms of an invariant associated with Equation 3 of [doi:10.1002/chem.200902300](http://dx.doi.org/10.1002/chem.200902300). The distinction between the general rhombic-eigenvalue case and the constructed axial, zero-trace case is explicit in the source.

## Setup and checks

Run `invariants()`. In each of 100 rhombic trials the function forms a random real symmetric 3×3 matrix as `A=randn(3); A=A+A'`, takes its eigenvalues and shuffles them. With `Ax=2*D(3)-(D(1)+D(2))` and `Rh=D(1)-D(2)`, it compares `(Ax^2+3*Rh^2)/6` against `(2/3)*(sum(D.^2)-D(1)*D(2)-D(1)*D(3)-D(2)*D(3))`; the absolute difference must be at most `1e-6`.

A separate 100-trial axial check starts from two random values with the third set equal to the second, subtracts the mean to make the triple zero-trace, and shuffles it. It compares `(D(1)-D(2))^2` with `sum(D.^2)-D(1)*D(2)-D(1)*D(3)-D(2)*D(3)`, again using an absolute `1e-6` tolerance. These are algebraic checks on sampled eigenvalue triples, not relaxation simulations.

## Observable result and scope

The first failed comparison raises its corresponding rhombic or axial test error. On completing both loops, the function prints `Axial test passed.`; it creates no plot and returns no fitted parameter. No random seed is set, so the tested triples vary with MATLAB's random stream.
