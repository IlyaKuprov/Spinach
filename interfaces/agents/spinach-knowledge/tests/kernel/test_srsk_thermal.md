# tests/kernel/test_srsk_thermal.m

- Source: `tests/kernel/test_srsk_thermal.m`
- Signature: `result=test_srsk_thermal()`
- Total lines: 167

## Purpose

Once-only thermalisation of additive SRSK relaxation and once-only addition of bosonic mode dissipation.

## Physical / mathematical content

Fast-source scalar relaxation rates and the thermal stationary state are checked with an oriented noncommuting quadrupolar interaction. A spectator cavity distinguishes unital dephasing from non-unital amplitude damping at finite temperature and checks trace conservation.

## Numerical / algorithmic content

Compares zero-destination, IME, and DiBari construction with explicit rate augmentation and once-only references across coupling signs and supported retention policies. Spin-boson cases cover damping, dephasing, both, and neither, at nonzero and zero scalar coupling, against spin-only thermalisation plus one original-temperature mode dissipator and the no-SRSK path with explicitly augmented rates.

## Syntax

`result=test_srsk_thermal()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
