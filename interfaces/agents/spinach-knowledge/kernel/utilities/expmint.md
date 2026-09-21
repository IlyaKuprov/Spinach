# kernel/utilities/expmint.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/expmint.m`
- Signature: `R=expmint(spin_system,A,B,C,T)`
- Total lines: 91

## Purpose

Computes matrix exponential integrals of the following general type: Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}] Matrix A must be Hermitian. For further info see the paper by Char- les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax: R=expmint(spin_system,A,B,C,T)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A,B,C -the three matrices involved in the integral
- T -integration time
- Output:
- R -the resulting integral
- Note: the auxiliary matrix method is massively faster than either
- commutator series or diagonalisation.
- Note: this is the most memory-intensive stage in a lot of calcula-
- tions; memory recycling is aggressive.

## Implementation structure

- Computes matrix exponential integrals of the following general type:
- Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]
- Matrix A must be Hermitian. For further info see the paper by Char-
- les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
- R=expmint(spin_system,A,B,C,T)
- A,B,C -the three matrices involved in the integral
- T -integration time
- Output:
- R -the resulting integral
- Note: the auxiliary matrix method is massively faster than either
- commutator series or diagonalisation.
- Note: this is the most memory-intensive stage in a lot of calcula-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nnz()`, `spalloc()`, `spdiags()`, `clear()`, `propagator()`, `report()`, `clean_up()`, `ismatrix()`, `all()`, `ishermitian()`, `isscalar()`.
