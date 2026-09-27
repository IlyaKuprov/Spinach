# examples/fundamentals/operator_tests/commutation_3.m

- Signature: `commutation_3()`

## Purpose

Checks single-spin and product-operator commutators across three Spinach formalisms.

## Physical / mathematical content

The model contains 1H and 235U spins at a 14.1 field setting, with scalar Zeeman entries 2.5 and 1.0 and symmetric scalar coupling entries of 10. In each of `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, it verifies the 235U ladder and Cartesian angular-momentum relations. It also constructs the two-spin products `L+`/`L+` and `L-`/`L-` and checks their commutators with the 235U and 1H `Lz` operators.

## Numerical / algorithmic content

Seven Frobenius-norm residuals are recorded for each formalism. The example passes when the Frobenius norm of the full residual array is below `1e-6`; otherwise it raises an error.

## Implementation structure

For each formalism, the script builds the spin system and basis, obtains the required single- and two-spin operators, evaluates the seven identities, and reports one cross-formalism result.
