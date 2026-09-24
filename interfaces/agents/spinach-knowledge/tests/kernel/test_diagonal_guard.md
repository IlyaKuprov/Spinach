# tests/kernel/test_diagonal_guard.m

- Source: `tests/kernel/test_diagonal_guard.m`
- Signature: `result=test_diagonal_guard()`
- Total lines: 103

## Purpose

Supported formalism boundaries for diagonal relaxation retention.

## Physical / mathematical content

Trace preservation, identity stationarity, Hermiticity, and longitudinal/transverse decay rates are checked for spin-half and spin-one baths.

## Numerical / algorithmic content

Requires explicit rejection of Zeeman diagonal retention while retaining spherical diagonal, Zeeman full, and omitted-retention controls.

## Syntax

`result=test_diagonal_guard()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
