# tests/kernel/test_cce_clusters.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_clusters.m)

## Purpose

Check frontier enumeration, closure, powers and physical group costs.

## Behaviour

Checks the following named assertions in the source:

- `Brute-force families through seven vertices`.
- `Exact Mobius incidence identities`.
- `Every labelled graph through four vertices`.
- `Explicit closure and cap overrun`.
- `Composite physical Hilbert dimension`.
- `Conservative central-mediated edge`.
- `Explicit approximate graph`.
- `Central-mediated full-Hilbert reference`.
- `Central mediation is not disconnected factorisation`.
- `Physical operator norm`.
- `Hamiltonian cutoff independence`.
- `Explicit resonance score`.
- `Union-of-centres bath radius`.
- `Deterministic score cap`.
- `Lexicographic score ties`.
- `One-group exact physical cost`.
- `Sparse large-bath capped enumeration`.
- `Duplicate explicit-family refusal`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_clusters()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_clusters.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
