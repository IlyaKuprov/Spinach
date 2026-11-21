# tests/kernel/test_mecce_exact.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_exact.m)

## Purpose

Zero jumps, dense Liouville full bath, and correlated support.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.
- `ME-gCCE Hahn zero jumps`.
- `ME-CCE Hahn zero jumps`.
- `Full-bath projected ME reference`.
- `Dissipative Hahn reference`.
- `C7 contract check`.
- `PyCCE reconstructed ME-CCE2`.
- `PyCCE reconstructed ME-gCCE2`.
- `ME tutorial validity`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_mecce_exact()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_exact.m)

Supplied-spin_config regression checks compare gCCE and ME-gCCE density, signal, and observables (raw and normalised) with an independent diagonal Hamiltonian solution. They also verify projected bath phases and swapped user projections, and require informative rejection of unsupported spin_config state/projection/reference aliases.

## Native contract

The registered test uses clexpn with free_evolution and hahn_echo, natural numeric density/signal outputs, and native compile/subsystem/generator/evaluator services. Its selective bath pulse uses refocusing.spins=2, operator=pi*Sx, and fraction=0.4. All independent dense Liouville, zero-jump, correlated-support, PyCCE-reference, normalisation, memory-budget, and configuration-rejection assertions remain. Observed preparation belongs to parameters.rho0; bath configurations contain bath densities only, so the state-override checks bind their explicit initial observed density before direct evaluation. Legacy initial_rho, projection, and reference aliases are rejected at the corresponding native boundary.

The independent configuration-state and reconstructed PyCCE tutorial fixtures execute inline with separate named physical and parameter structures, preserving the surrounding exact-dynamics fixture.
