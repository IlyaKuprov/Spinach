# tests/kernel/test_keyhole_guard.m

- Source: `tests/kernel/test_keyhole_guard.m`
- Signature: `result=test_keyhole_guard()`
- Total lines: 163

## Purpose

Explicit boundaries for keyhole-projected exact-Hessian optimisation.

## Physical / mathematical content

Noncommuting spin-half controls and population projections test supported gradients and Hessians against finite differences.

## Numerical / algorithmic content

Checks setup and direct-engine refusals, retained method identities, first-order keyholes, empty schedules, and supported Hilbert counterparts.

## Syntax

`result=test_keyhole_guard()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
