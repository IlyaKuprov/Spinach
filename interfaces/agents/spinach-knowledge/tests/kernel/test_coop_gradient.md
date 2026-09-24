# tests/kernel/test_coop_gradient.m

- Source: `tests/kernel/test_coop_gradient.m`
- Signature: `result=test_coop_gradient()`
- Total lines: 219

## Purpose

Cooperative phase gradients for the requested transfer fidelity.

## Physical / mathematical content

The cooperative objective combines primary transfer and squared orthogonal impurity; independent density-matrix propagation reconstructs that objective.

## Numerical / algorithmic content

Finite differences check phase derivatives in both density formalisms with unit/nonunit complex targets, noncommuting pulses, and a power ensemble.

## Syntax

`result=test_coop_gradient()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
