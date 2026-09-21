# kernel/utilities/rspert.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rspert.m`
- Signature: `[Ep,Vp]=rspert(E0,H1,order)`
- Total lines: 105

## Purpose

Rayleigh-Schrodinger perturbation theory to arbitrary order, Eqs 2.21-2.23 from Stefan Stoll's PhD thesis, with the typo fixed in the numerator of Eq 2.21. Syntax: [Ep,Vp]=rspert(E0,H1,order)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used, 6
- is the sensible maximum

## Outputs

- Ep -eigenvalues of H0+H1 to the specified order,
- a vector of reals, not necessarily sorted in
- the same way as the input
- Vp -normalised eigenvectors of H0+H1 to the spe-
- cified order in perturbation theory, a squa-
- re unitary matrix with eigenvectors in cols
- in the same order as the eigenvalues in Ep
- Notes: there must be no degeneracies in H0; H1 must be Hermitian,
- the theory only converges when norm(H1,2) is much smaller
- than the smallest energy gap in H0; numerical artefacts
- appear beyond sixth order; complexity is linear in the or-
- der and cubic in the matrix dimension.

## Implementation structure

- Rayleigh-Schrodinger perturbation theory to arbitrary order, Eqs
- 2.21-2.23 from Stefan Stoll's PhD thesis, with the typo fixed in
- the numerator of Eq 2.21. Syntax:
- [Ep,Vp]=rspert(E0,H1,order)
- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used, 6
- is the sensible maximum
- Ep -eigenvalues of H0+H1 to the specified order,
- a vector of reals, not necessarily sorted in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `logical()`, `speye()`, `any()`, `iscolumn()`, `ishermitian()`, `isscalar()`.
