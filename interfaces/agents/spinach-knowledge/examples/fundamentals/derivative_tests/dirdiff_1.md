# examples/fundamentals/derivative_tests/dirdiff_1.m

- Signature: `dirdiff_1()`

## Purpose

Check first- and second-order derivatives of the matrix exponential against finite differences.

## Physical / mathematical content

The test exercises derivatives of a propagator with respect to one or two Hermitian matrix directions. It loops over three Spinach formalisms to construct test systems; the derivative comparison itself uses random matrices and is formalism-independent.

## Numerical / algorithmic content

A centered two-point difference checks the first derivative with step `1e-3`; a centered four-point mixed difference checks the second derivative. The numerical and analytical second derivatives are compared after averaging the two perturbation orders. Relative-error limits are `1e-5` (first derivative) and `1e-3` (second derivative).

## Implementation structure

- Construct a test Spinach system for each of `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`.
- Generate a random `5×5` Hermitian Hamiltonian and two Hermitian direction matrices.
- Compare finite differences of `propagator(...,1)` with the first- and second-derivative results from `dirdiff`.
- Raise an error if either relative-error limit is exceeded.
