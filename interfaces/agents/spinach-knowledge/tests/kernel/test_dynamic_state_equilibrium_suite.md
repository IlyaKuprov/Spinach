# tests/kernel/test_dynamic_state_equilibrium_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_state_equilibrium_suite.m`
- Signature: `result=test_dynamic_state_equilibrium_suite()`
- Total lines: 73

## Purpose

Tests thermal equilibrium state construction paths. Syntax: result=test_dynamic_state_equilibrium_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- result -regression test result with explanatory messages
- The test checks Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian
- thermal equilibrium construction against direct Boltzmann references.

## Implementation structure

- Tests thermal equilibrium state construction paths. Syntax:
- result=test_dynamic_state_equilibrium_suite()
- result -regression test result with explanatory messages
- The test checks Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian
- thermal equilibrium construction against direct Boltzmann references.
- Announce the test target
- State the equilibrium-constructor target of the test
- Build a one-spin Hilbert-space system with finite temperature
- Set an explicit non-degenerate Hamiltonian in angular frequency units
- Compare Hilbert-space equilibrium to the direct Boltzmann density matrix
- Build the matching Zeeman-Liouville system
- Compare left-product Liouville equilibrium to the vectorised density matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `equilibrium()`, `test_spin_system()`, `test_close()`, `speye()`, `rho_ref()`, `deal()`, `orientation()`.
