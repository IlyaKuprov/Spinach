# tests/kernel/test_giant_ham_descr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_giant_ham_descr.m`
- Signature: `result=test_giant_ham_descr()`
- Total lines: 134

## Purpose

Tests the giant spin Hamiltonian descriptor route. Syntax: result=test_giant_ham_descr()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `local_case()`, `local_ref()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks that high-rank giant spin Hamiltonian terms assembled
- through the descriptor route match direct spherical-tensor assembly.

## Implementation structure

- Tests the giant spin Hamiltonian descriptor route. Syntax:
- result=test_giant_ham_descr()
- result -regression test result with explanatory messages
- The test checks that high-rank giant spin Hamiltonian terms assembled
- through the descriptor route match direct spherical-tensor assembly.
- Announce the test target
- State the Hamiltonian target of the test
- Check complete giant spin terms
- Check secular giant spin terms
- Checks one giant spin Hamiltonian assumption
- Build a compact high-rank giant spin system
- Apply the requested Hamiltonian assumption

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_case()`, `test_spin_system()`, `assume()`, `hamiltonian()`, `orientation()`, `local_ref()`, `test_close()`, `double()`, `mprealloc()`, `wigner()`, `euler_angles()`, `num2str()`, `operator()`.
