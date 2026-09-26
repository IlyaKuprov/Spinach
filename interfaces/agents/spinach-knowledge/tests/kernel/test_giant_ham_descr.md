# tests/kernel/test_giant_ham_descr.m

- Signature: `result=test_giant_ham_descr()`

## Purpose

Tests the giant spin Hamiltonian descriptor route. Syntax: result=test_giant_ham_descr()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
