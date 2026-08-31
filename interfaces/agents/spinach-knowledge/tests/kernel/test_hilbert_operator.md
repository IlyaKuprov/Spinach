# tests/kernel/test_hilbert_operator.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_hilbert_operator.m`
- Signature: `result=test_hilbert_operator()`
- Total lines: 48

## Purpose

Tests Hilbert-space operator generation. Syntax: result=test_hilbert_operator()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Hilbert-space operator generation\n')`.
- Lines 19-22: State the physical target of the test; implemented by `result=new_test_result('kernel/hilbert_operator', 'Hilbert-space operator generation', 'operator() must map Lx, Ly, Lz, L+, and L- labels to spin matrices.')`.
- Lines 24-25: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 32-33: Textbook spin-half reference matrices; implemented by `S=pauli(2)`.
- Lines 35-37: Check label-to-matrix mapping; implemented by `result=test_close(result,'Lx operator',operator(spin_system,'Lx',1),S.x,1e-15,1e-15, 'a one-proton Lx operator is the spin-half Sx matrix')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/hilbert_operator', 'Hilbert-space operator generation', 'operator() must map Lx, Ly, Lz, L+, and L- labels to spin matrices.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `S` using `S=pauli(2)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that operator() builds the correct one-spin Hilbert-space
- angular momentum matrices from human-readable labels.

## Implementation structure

- Tests Hilbert-space operator generation. Syntax:
- result=test_hilbert_operator()
- result -regression test result with explanatory messages
- The test checks that operator() builds the correct one-spin Hilbert-space
- angular momentum matrices from human-readable labels.
- Announce the test target
- State the physical target of the test
- Build a one-proton Hilbert-space spin system
- Textbook spin-half reference matrices
- Check label-to-matrix mapping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `operator()`, `test_spin_system()`, `pauli()`, `test_close()`.
