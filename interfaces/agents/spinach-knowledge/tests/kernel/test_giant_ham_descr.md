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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Giant spin Hamiltonian descriptor\n')`.
- Lines 19-22: State the Hamiltonian target of the test; implemented by `result=new_test_result('kernel/giant_ham_descr', 'Giant spin Hamiltonian descriptor', 'high-rank giant spin terms must match direct spherical-tensor assembly.')`.
- Lines 24-25: Check complete giant spin terms; implemented by `result=local_case(result,'labframe','strong',[0.41 0.29 0.13])`.
- Lines 27-28: Check secular giant spin terms; implemented by `result=local_case(result,'deer-zz','secular',[0.17 0.39 0.51])`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/giant_ham_descr', 'Giant spin Hamiltonian descriptor', 'high-rank giant spin terms must match direct spherical-tensor assembly.')`.

### Local helper functions

- Line 34: `local_case()` — `function result=local_case(result,assumption,strength,euler_angles)`. Build a compact high-rank giant spin system
  - Representative operation: `sys.magnet=0`.
  - Representative operation: `sys.isotopes={'E8'}`.
- Line 74: `local_ref()` — `function H_ref=local_ref(spin_system,euler_angles)`. Start from a zero Hamiltonian
  - Representative operation: `H_ref=mprealloc(spin_system,0)`.
  - Representative operation: `for n=1:spin_system.comp.nspins`.

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
