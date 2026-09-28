# tests/kernel/test_homodec_adapt.m

- Signature: `result=test_homodec_adapt()`

## Purpose

Acquisition proton irradiation through Hilbert formalism admission.

## Physical / mathematical content

A phase-shifted soft pulse and noncommuting proton irradiation in a 1H-13C system must agree with native Liouville dynamics.

## Numerical / algorithmic content

Checks operator conversion, nonzero signals and irradiation effects, zero/absent power fields, dead time, and no-op formalisms.

## Syntax

`result=test_homodec_adapt()`

## Parameters / inputs

None. The test constructs its own physical system.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

Optional operator caching is not enabled.
