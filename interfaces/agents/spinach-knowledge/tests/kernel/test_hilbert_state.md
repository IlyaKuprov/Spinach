# tests/kernel/test_hilbert_state.m

- Signature: `result=test_hilbert_state()`

## Purpose

Tests that `state()` returns the expected Hilbert-space density matrices for a one-spin system.

## Test setup

Creates a one-proton (`1H`) spin system with zero magnetic field (`sys.magnet=0`), zero scalar Zeeman interaction, `zeeman-hilb` formalism, and no basis approximation. Reference spin-half matrices are obtained from `S=pauli(2)`.

## Assertions

Each comparison uses `test_close` with absolute and relative tolerances of `1e-15`:

- `state(spin_system,'Lz',1)` equals `S.z` (longitudinal magnetisation).
- `state(spin_system,'Lx',1)` equals `S.x` (transverse in-phase density matrix).
- `state(spin_system,'E',1)` equals `S.u` (unit density matrix).

## Output

`result` is a regression test result with explanatory messages, created for `kernel/hilbert_state` and updated by the three comparisons.