# examples/fundamentals/derivative_tests/dirdiff_2.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_2.m)

- Signature: `dirdiff_2()`

## Purpose

This example compares the analytical left- and right-control derivatives returned by `trapdiff` for the second-order Magnus product quadrature with centred finite differences of the associated matrix exponential.

## Test construction

For each of `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`, the script constructs a test spin system; the derivative comparison itself uses independent random complex 50-by-50 drift matrices and a 50-by-50 control matrix, not a physical spin Hamiltonian. It scales the time step with the reciprocal drift-matrix norms and tests the coherent plus non-symmetric dissipative case.

## Derivatives and comparison

It uses a step of `sqrt(eps('double'))`, constructs left and right directions for the second-order Magnus trapezoid, and estimates each matrix-exponential derivative by central differences in that direction. It compares both estimates with `trapdiff` using spectral 2-norm residuals; each must be strictly below `10*sqrt(eps('double'))`. Otherwise the script raises an error; a successful comparison prints `trapezium quadrature derivative test passed`.

## Scope

This is a finite-difference consistency check for the stated matrix construction and the two derivative outputs. The formalism labels do not constitute three independent derivative calculations, and the random matrices do not by themselves establish behaviour for every physical generator or for a composed simulation algorithm.