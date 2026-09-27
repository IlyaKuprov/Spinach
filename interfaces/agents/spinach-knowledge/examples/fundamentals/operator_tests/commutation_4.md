# examples/fundamentals/operator_tests/commutation_4.m

- Signature: `commutation_4()`

## Purpose

Checks the central-transition operator algebra in three Spinach representations.

## Physical / mathematical content

The system specification includes 1H and 235U spins, scalar Zeeman entries 2.5 and 1.0, and symmetric scalar coupling entries of 10. For each of `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, the script constructs the 235U central-transition Cartesian and ladder operators and tests `[CTz,CT+]=CT+`, `[CTz,CT-]=-CT-`, and `[CTx,CTy]=iCTz`.

## Numerical / algorithmic content

The three Frobenius-norm residuals per formalism are collected in an array. The test passes if its total Frobenius norm is below `1e-6`; otherwise it raises an error.

## Implementation structure

The spin system and basis are rebuilt for each formalism, the five central-transition operators are obtained, and the three commutator identities are checked before a combined pass/failure report.
