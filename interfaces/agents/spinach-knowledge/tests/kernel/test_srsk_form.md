# tests/kernel/test_srsk_form.m

- Signature: `result=test_srsk_form()`

## Purpose

SRSK formalism refusal and supported fast-source relaxation.

## Physical / mathematical content

A rapidly relaxing nitrogen source produces the expected proton scalar-relaxation rates in spherical-tensor Liouville space.

## Numerical / algorithmic content

Checks explicit unsupported Zeeman requests across retention/equilibrium choices and preserves Zeeman Lindblad and spherical SRSK controls.

## Syntax

`result=test_srsk_form()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
