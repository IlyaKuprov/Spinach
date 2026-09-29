# tests/kernel/test_hilbert_operator.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilbert_operator.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilbert_operator.m)

## Purpose

Regression test for Hilbert-space operator generation. It verifies that `operator()` builds the correct one-spin Hilbert-space angular momentum matrices from human-readable labels.

## Behaviour

- Announces the test target with `fprintf('TESTING: Hilbert-space operator generation\n')`.
- Initialises a test result via `new_test_result('kernel/hilbert_operator', 'Hilbert-space operator generation', 'operator() must map Lx, Ly, Lz, L+, and L- labels to spin matrices.')`.
- Builds a one-proton Hilbert-space spin system with:
  - `sys.magnet = 0`
  - `sys.isotopes = {'1H'}`
  - `inter.zeeman.scalar = {0}`
  - `bas.formalism = 'zeeman-hilb'`
  - `bas.approximation = 'none'`
  - The spin system is created with `test_spin_system(sys, inter, bas)`.
- Uses `pauli(2)` to obtain textbook spin-half reference matrices.
- Checks label-to-matrix mapping with `test_close` using tolerances `1e-15` (absolute) and `1e-15` (relative):
  - `operator(spin_system, 'Lx', 1)` against `S.x` — 'a one-proton Lx operator is the spin-half Sx matrix'
  - `operator(spin_system, 'Ly', 1)` against `S.y` — 'a one-proton Ly operator is the spin-half Sy matrix'
  - `operator(spin_system, 'Lz', 1)` against `S.z` — 'a one-proton Lz operator is the spin-half Sz matrix'
  - `operator(spin_system, 'L+', 1)` against `S.p` — 'L+ is the spin raising operator in the Zeeman basis'
  - `operator(spin_system, 'L-', 1)` against `S.m` — 'L- is the spin lowering operator in the Zeeman basis'

## Inputs and outputs

- **Syntax:** `result = test_hilbert_operator()`
- **Outputs:**
  - `result` — regression test result with explanatory messages.
- **Inputs:** None.

## References

- Source file: [tests/kernel/test_hilbert_operator.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_hilbert_operator.m) on GitHub.
