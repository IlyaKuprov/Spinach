# tests/kernel/test_hilbert_state.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_hilbert_state.m`
- Signature: `result=test_hilbert_state()`
- Total lines: 44

## Purpose

Tests Hilbert-space state generation. Syntax: result=test_hilbert_state()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Hilbert-space state generation\n')`.
- Lines 19-22: State the physical target of the test; implemented by `result=new_test_result('kernel/hilbert_state', 'Hilbert-space state generation', 'state() must map observable labels to density matrices.')`.
- Lines 24-25: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 32-33: Textbook spin-half reference matrices; implemented by `S=pauli(2)`.
- Lines 35-37: Check density matrices generated from state labels; implemented by `result=test_close(result,'Lz state',state(spin_system,'Lz',1),S.z,1e-15,1e-15, 'the Lz state is the longitudinal magnetisation density matrix')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/hilbert_state', 'Hilbert-space state generation', 'state() must map observable labels to density matrices.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `S` using `S=pauli(2)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that state() returns the expected density matrices in
- Hilbert space for a one-spin system.

## Implementation structure

- Tests Hilbert-space state generation. Syntax:
- result=test_hilbert_state()
- result -regression test result with explanatory messages
- The test checks that state() returns the expected density matrices in
- Hilbert space for a one-spin system.
- Announce the test target
- State the physical target of the test
- Build a one-proton Hilbert-space spin system
- Textbook spin-half reference matrices
- Check density matrices generated from state labels

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `state()`, `test_spin_system()`, `pauli()`, `test_close()`.
