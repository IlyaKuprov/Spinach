# tests/kernel/test_dynamic_state_equilibrium_suite.m

- Signature: `result=test_dynamic_state_equilibrium_suite()`

## Purpose

Regression tests for equilibrium-state construction.

## Tests

- Compares `equilibrium()` with the explicitly normalized Boltzmann matrix `expm(-beta*H)/trace(expm(-beta*H))`.
- Checks the Liouville-space state representation.
- Compares the Euler-angle Hamiltonian route `equilibrium(spin_h,H,Q,euler_angles)` with the equivalent call using the pre-oriented Hamiltonian, `equilibrium(spin_h,H_oriented)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian
- thermal equilibrium construction against direct Boltzmann references.
