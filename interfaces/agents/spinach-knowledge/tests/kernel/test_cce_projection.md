# tests/kernel/test_cce_projection.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_projection.m)

## Purpose

Check explicit states, projected fields, and resolved virtual levels.

## Behaviour

Checks the following named assertions in the source:

- `alpha expectation`.
- `beta expectation`.
- `second order off by default`.
- `virtual coefficient`.
- `second-order spectrum`.
- `near-degenerate denominator rejected`.
- `overlap guard`.
- `gCCE correction rejected`.
- `multiple observed-spin system`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_projection()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_projection.m)

Checks second-order mediated coupling-graph dynamics against independently constructed dense projected Hamiltonians, and verifies that omitting induced edges requires the explicit graph approximation flag.

The mediated-pair and explicit graph-approximation cases use `clexpn(system,@free_evolution,parameters,'labframe')`. The pair uses the explicit ordered magnetic quantum numbers `[0.5 -0.5]`, which select the same two Zeeman basis columns as the direct projection checks. Physical observed IDs and transition selectors are experiment parameters; second-order, graph, bath, order, and numerical controls belong to `parameters.cce`. The native `cce_compile` service supplies both physical models used by direct internal `cce_project` checks. The internal projection control structure is distinct from the context parameters; it preserves explicit-state, gap, overlap, and double-counting rejection assertions.
