# tests/kernel/test_dynamic_state_projection_suite.m

- Signature: `result=test_dynamic_state_projection_suite()`

## Purpose

Regression tests for spin-state projections and related operator/state diagnostics.

## Tests

- Checks coherence operators and their adjoints for selected coherence orders and spins.
- Checks commutators involving dephased populations.
- Compares captured `stateinfo()` text with the expected report.
- Tests projection onto an isotropic zero-field triplet state.

## Outputs

- result -regression test result with explanatory messages
- The test checks deuteron-pair coherences, dephased population stationarity,
- captured stateinfo() output, and isotropic zero-field triplet projection.
