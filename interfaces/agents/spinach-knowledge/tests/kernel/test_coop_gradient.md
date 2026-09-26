# tests/kernel/test_coop_gradient.m

- Source: `tests/kernel/test_coop_gradient.m`
- Signature: `result=test_coop_gradient()`
- Total lines: 318

## Purpose

Cooperative phase gradients for the requested transfer fidelity.

## Physical / mathematical content

The cooperative objective combines primary transfer and squared orthogonal impurity; independent density-matrix propagation reconstructs that objective.

## Numerical / algorithmic content

Finite differences check phase derivatives in both density formalisms with unit/nonunit complex targets, noncommuting pulses, a power ensemble, vanishing impurity, and purely imaginary auxiliary overlaps. Exact zero auxiliary values and gradients remain valid engine outputs. A cooperative one-spin transfer with nonzero primary overlap and gradient but near-zero combined score must be admitted and improve in a real optimiser step.

All four `fmaxnewton` methods are tested in both density formalisms: unusable assembled initial guesses and all-frozen gradients must receive the poor-guess diagnostic without singular-solve warnings. Constant objectives remain valid for zero-iteration evaluation. A nonstationary physical transfer must improve over multiple iterations, with Hessian evaluations requested only by Newton and Goodwin, once per iteration.

## Syntax

`result=test_coop_gradient()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest`; its tiny optimiser checks are correctness tests, not performance benchmarks.
