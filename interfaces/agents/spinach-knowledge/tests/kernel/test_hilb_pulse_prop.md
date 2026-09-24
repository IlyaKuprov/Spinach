# tests/kernel/test_hilb_pulse_prop.m

- Source: `tests/kernel/test_hilb_pulse_prop.m`
- Signature: `result=test_hilb_pulse_prop()`
- Total lines: 246

## Purpose

Regresses the reusable one-sided propagator returned by `shaped_pulse_xy` in Hilbert space while preserving its two-sided density-state evolution.

## Physical / mathematical content

For every method and quadrature, the returned operator is the ordered product of slice exponentials, with later slices multiplying on the left. It must satisfy `rho_final=P*rho_initial*P'`, preserve unitarity, density trace, and Hermiticity, and work for another initial density. Complex noncommuting spin-half generators expose ordering errors; an identity contribution to the Hamiltonian also tests the propagator phase that cancels from density evolution.

## Numerical / algorithmic content

All six `expv-*`, `expm-*`, and `evol-*` choices are compared with independently assembled dense exponentials and two-point Lie generators. Nonconstant pulses, constant pulses, zero durations, both small-matrix and commutator-series Hilbert state paths, every trajectory point, and output-count invariance are covered. Wavefunction and Liouville calls provide one-sided controls.

Dimension-512 direct-sum fixtures run all six choices with sparse and full storage on the CPU and with full storage on an available GPU. Sparse GPU coverage includes `expm-pwc` and `evol-pwc` unconditionally; the other four methods require MATLAB sparse-GPU scalar division, whose availability is probed directly before their execution. An unavailable operation produces an explicit per-method `SKIP`, not a claimed pass. The matrix-density threshold keeps the two storage routes distinct; the dimension exceeds the explicit propagator GPU dispatch threshold. Analytic block references check the operator, state, trajectory, invariants, and final host-memory outputs. GPU checks emit an explicit `SKIP` message when `canUseGPU` is unavailable or false; errors in an available GPU production path are not caught or converted into skips.

## Syntax

```matlab
result=test_hilb_pulse_prop()
```

## Parameters / inputs

None. Requires the Spinach production and regression-test libraries on the MATLAB path. GPU coverage additionally requires Parallel Computing Toolbox and a usable GPU.

## Outputs

`result` contains the regression check messages and failures. A CPU pass with a GPU skip is not evidence that the GPU path passed.

## Header notes

The test is registered as `kernel/hilb_pulse_prop` in `tests/lib/test_manifest.m`.
