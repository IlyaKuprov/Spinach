# tests/kernel/test_hilbert_state.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilbert_state.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilbert_state.m)

## Purpose

Regression test for Hilbert-space state generation. It verifies that `state()` maps observable labels to the expected density matrices for a one-spin system.

## Behaviour

- Announces the test target with `fprintf('TESTING: Hilbert-space state generation\n')`.
- Initialises a test result via `new_test_result('kernel/hilbert_state', 'Hilbert-space state generation', 'state() must map observable labels to density matrices.')`.
- Builds a one-proton Hilbert-space spin system with:
  - `sys.magnet = 0`
  - `sys.isotopes = {'1H'}`
  - `inter.zeeman.scalar = {0}`
  - `bas.formalism = 'zeeman-hilb'`
  - `bas.approximation = 'none'`
  - the system is produced by `test_spin_system(sys, inter, bas)`.
- Obtains textbook spin-half reference matrices from `pauli(2)`.
- Runs three closeness checks with `test_close`, each using absolute and relative tolerances of `1e-15`:
  - `'Lz state'`: `state(spin_system,'Lz',1)` against `S.z`, described as the longitudinal magnetisation density matrix.
  - `'Lx state'`: `state(spin_system,'Lx',1)` against `S.x`, described as the transverse in-phase density matrix.
  - `'identity state'`: `state(spin_system,'E',1)` against `S.u`, described as the unit density matrix.

## Inputs and outputs

- **Outputs:**
  - `result` — regression test result with explanatory messages.
- **Inputs:** none; the function is called as `result = test_hilbert_state()`.

## References

- Source file header documents the syntax `result=test_hilbert_state()` and the output description.
