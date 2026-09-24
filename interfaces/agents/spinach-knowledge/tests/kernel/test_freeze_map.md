# tests/kernel/test_freeze_map.m

- Source: `tests/kernel/test_freeze_map.m`
- Signature: `result=test_freeze_map()`
- Total lines: 157

## Purpose

Frozen input derivatives through composed waveform transformations.

## Physical / mathematical content

The chain rule must include every physical waveform coordinate before constraining the input coordinates.

## Numerical / algorithmic content

Checks serial filters, phase rotation, power scaling, empty masks, supported exact Hessians, dissipative dynamics, and unchanged direct-engine masking.

## Syntax

`result=test_freeze_map()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
