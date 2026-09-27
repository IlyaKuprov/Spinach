# tests/kernel/test_homodec_adapt.m

- Signature: `result=test_homodec_adapt()`

## Purpose

Acquisition irradiation through Hilbert formalism admission.

## Physical / mathematical content

A phase-shifted soft pulse and noncommuting heteronuclear irradiation must agree with native Liouville dynamics.

## Numerical / algorithmic content

Checks operator conversion, nonzero signals and irradiation effects, zero/absent power fields, dead time, and no-op formalisms.

## Syntax

`result=test_homodec_adapt()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
