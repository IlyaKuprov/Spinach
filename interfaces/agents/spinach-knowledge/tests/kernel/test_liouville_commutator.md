# tests/kernel/test_liouville_commutator.m

- Signature: `result=test_liouville_commutator()`

## Purpose

Tests Liouville-space commutation superoperators by checking that `operator()` in `'comm'` mode equals its `'left'` mode minus its `'right'` mode.

## Physical / mathematical content

- For an operator `A` acting on a density operator `rho`, the commutation superoperator is `A*rho-rho*A`.

## Numerical / algorithmic content

- Constructs a one-proton spin system in the `'zeeman-liouv'` formalism with no basis approximation, zero magnetic field, and zero scalar Zeeman interaction.
- Generates the `Lz` superoperator for spin 1 in `'comm'`, `'left'`, and `'right'` modes, then checks `L_comm = L_left-L_right` with absolute and relative tolerances of `1e-15`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the test and initializes its result record.
- Builds the spin system, generates the three superoperators, and records the identity check using `test_close()`.
