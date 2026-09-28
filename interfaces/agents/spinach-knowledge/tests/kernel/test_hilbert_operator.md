# tests/kernel/test_hilbert_operator.m

- Signature: `result=test_hilbert_operator()`

## Purpose

Tests that `operator()` maps human-readable labels to the correct one-spin Hilbert-space angular momentum matrices.

## Test setup

- Builds a one-proton (`1H`) spin system with zero magnetic field and zero scalar Zeeman interaction.
- Uses the `zeeman-hilb` formalism with no approximation.
- Obtains spin-half reference matrices from `pauli(2)`.

## Checks

Each `operator()` result is compared with its reference matrix using `test_close` with both tolerances set to `1e-15`:

| Label | Reference matrix |
| --- | --- |
| `Lx` | `S.x` |
| `Ly` | `S.y` |
| `Lz` | `S.z` |
| `L+` | `S.p` (raising operator) |
| `L-` | `S.m` (lowering operator) |

## Output

- `result` — regression test result with explanatory messages.