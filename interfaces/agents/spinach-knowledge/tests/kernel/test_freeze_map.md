# tests/kernel/test_freeze_map.m

Regression test for frozen input derivatives in Spinach optimal-control waveform maps.

Source: [tests/kernel/test_freeze_map.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_freeze_map.m)

## Purpose

`test_freeze_map` verifies that freezing waveform coordinates during optimal-control optimisation zeroes exactly the frozen derivatives while leaving the objective value and all unfrozen derivatives unchanged. The test exercises temporal filters, phase cycles, power scaling, and constrained Hessians in both the Hilbert and vectorised Liouville density formalisms, using noncommuting one-spin transfers.

## What is tested

The `control.freeze` mask applies to **input** waveform coordinates. Freezing an input must leave the objective value unchanged, set exactly its gradient component to zero, and preserve all free derivatives after pulling back through the physical waveform map. The suite checks this in Hilbert and vectorised Liouville formalisms with noncommuting spin transfers, including empty masks, temporal distortion, phase cycling, trapezium integration, phase-only controls, and coupled, rank-reducing or dimension-increasing coordinate maps. Free-gradient pullback is checked at `1e-12`; active derivatives are also compared with finite differences at `1e-8`.

For constrained curvilinear controls the same masking property covers every objective channel, including the norm-squared penalty. With exact-Hessian methods, frozen rows and columns of the Hessian must vanish and the active block must agree with finite differences of production gradients (`1e-8` tolerance). The direct `grape_liouv` mask remains effective; the direct `grape_hilb` engine retains its existing unmasked derivatives.

Finally, a long physically frozen interval must not require an unreachable propagation derivative. This remains true when a rank-reducing distortion and phase cycle together cancel the otherwise free input direction: the steady-state objective and gradient must stay finite, including with a non-Hermitian drift.

## Inputs and outputs

```matlab
result=test_freeze_map()
```

**Inputs:** none.

**Outputs:**

- `result` — regression test result structure with explanatory messages, accumulated through `test_close` and `test_true` assertions.

## References

- Spinach optimal control module functions used here: `optimcon`, `grape_xy`, `grape_curv`, `grape_phase`, `grape_liouv`, `grape_hilb`, `new_test_result`, `test_close`, `test_true`, `firf`, `spf`, `pauli`.
