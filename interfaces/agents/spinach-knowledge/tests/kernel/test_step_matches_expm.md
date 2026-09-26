# tests/kernel/test_step_matches_expm.m

- Signature: `result=test_step_matches_expm()`

## Purpose

Tests Hilbert-space propagation against matrix exponentiation. Syntax: result=test_step_matches_expm()

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- result -regression test result with explanatory messages
- The test checks the Spinach sign convention for density-matrix evolution:
- rho(t)=exp(-iHt) rho(0) exp(+iHt).

## Implementation structure

- Tests Hilbert-space propagation against matrix exponentiation. Syntax:
- result=test_step_matches_expm()
- result -regression test result with explanatory messages
- The test checks the Spinach sign convention for density-matrix evolution:
- rho(t)=exp(-iHt) rho(0) exp(+iHt).
- Announce the test target
- State the propagation target of the test
- Build a one-proton Hilbert-space spin system
- Define a Hamiltonian and an initial density matrix
- Build the independent exact propagator
- Check exact finite-dimensional propagation
