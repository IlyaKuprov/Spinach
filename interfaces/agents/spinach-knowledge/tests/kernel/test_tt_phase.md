# tests/kernel/test_tt_phase.m

- Source: `tests/kernel/test_tt_phase.m`
- Signature: `result=test_tt_phase()`
- Total lines: 98

## Purpose

Tensor-train absolute compression budgets under coefficient phases.

## Physical / mathematical content

A signed Hermitian spin interaction and genuinely complex cores test phase-independent Frobenius error bounds.

## Numerical / algorithmic content

Exercises multiple ranks/core counts, positive/negative/complex phases, actual rank reduction, zero tolerance, zero coefficients, and identity multiplication.

## Syntax

`result=test_tt_phase()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
