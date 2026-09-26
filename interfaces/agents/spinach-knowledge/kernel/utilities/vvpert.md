# kernel/utilities/vvpert.m

- Signature: `[Ep,G]=vvpert(E0,H1,order)`

## Purpose

Van Vleck perturbation theory, following Shavitt and Redmon, but excluding the quasi-degenerate split. Syntax: [Ep,G]=vvpert(E0,H1,order)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used; numerical
- artefacts appear beyond about 10-12 for typical
- problems

## Outputs

- Ep -eigenvalues of H0+H1 to the specified order,
- a column vector of reals, not necessarily
- sorted in the same way as the input
- G -Van Vleck transformation generator, such that
- expm(G) is a square unitary matrix with eigen-
- vectors in columns, in the same order as the
- eigenvalues in Ep
- Notes: there must be no degeneracies in H0; H1 must be Hermitian,
- the theory only converges when norm(H1,2) is much smaller
- than the smallest energy gap in H0; complexity is cubic
- both in the order and in the matrix dimension.

## Implementation structure

- Van Vleck perturbation theory, following Shavitt and Redmon, but
- excluding the quasi-degenerate split. Syntax:
- [Ep,G]=vvpert(E0,H1,order)
- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used; numerical
- artefacts appear beyond about 10-12 for typical
- problems
- Ep -eigenvalues of H0+H1 to the specified order,
- a column vector of reals, not necessarily
