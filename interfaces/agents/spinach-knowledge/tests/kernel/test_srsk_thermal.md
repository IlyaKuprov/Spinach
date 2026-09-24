# tests/kernel/test_srsk_thermal.m

- Source: `tests/kernel/test_srsk_thermal.m`
- Signature: `result=test_srsk_thermal()`
- Total lines: 106

## Purpose

Once-only thermalisation of additive SRSK relaxation.

## Physical / mathematical content

Fast-source scalar relaxation rates and the thermal stationary state are checked with an oriented noncommuting quadrupolar interaction.

## Numerical / algorithmic content

Compares zero-destination, IME, and DiBari construction with explicit rate augmentation and once-only references across coupling signs and supported retention policies.

## Syntax

`result=test_srsk_thermal()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
