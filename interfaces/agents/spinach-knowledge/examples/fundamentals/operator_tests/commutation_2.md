# examples/fundamentals/operator_tests/commutation_2.m

- Signature: `commutation_2()`

## Purpose

Checks angular-momentum commutation relations in three Spinach representations.

## Physical / mathematical content

The test builds a single 235U spin at zero field with zero scalar Zeeman interaction. For each of `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, it constructs the Cartesian and ladder operators and evaluates the three relations `[Lz,L+]=L+`, `[Lz,L-]=-L-`, and `[Lx,Ly]=iLz`.

## Numerical / algorithmic content

Each commutator residual is measured with the Frobenius norm. The three formalisms pass together only when the Frobenius norm of the 3-by-3 residual array is below `1e-6`; otherwise the example raises an error.

## Implementation structure

The script recreates the spin system and basis for each formalism, fills one residual column with the three operator identities, then reports a single cross-formalism pass or failure.
