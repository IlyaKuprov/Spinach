# tests/kernel/test_freeze_map.m

- Source: `tests/kernel/test_freeze_map.m`
- Signature: `result=test_freeze_map()`
- Total lines: 254

## Purpose

Frozen input derivatives through composed waveform transformations.

## Physical / mathematical content

The chain rule must include every physical waveform coordinate before constraining the input coordinates.

## Numerical / algorithmic content

Checks serial filters, phase rotation, power scaling, empty masks, supported exact Hessians, dissipative dynamics, and unchanged direct-engine masking.

The same two waveform/target fixtures are tested in Hilbert and vectorised Liouville formalisms with identity, coupled linear, three-to-two, one-to-two, and nonlinear polar coordinate maps. Empty and all-false masks, trapezium propagation, phase/power/filter composition, and a non-zero rectilinear norm-square penalty are covered. Frozen entries must vanish in every returned curvilinear gradient channel; free entries must equal the unmasked pullback and centred objective differences at two increments. Error and reference-derivative norms are reported. Phase-only gradient/Hessian masking and the existing direct Hilbert engine behaviour are checked in addition to the original Cartesian and direct Liouville cases.

## Syntax

`result=test_freeze_map()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
