# tests/kernel/test_giant_cache.m

- Signature: `result=test_giant_cache()`

## Purpose

Giant-spin Hamiltonian cache identity under retained interactions.

## Physical / mathematical content

Rotated Hermitian giant-spin Hamiltonians contain genuinely complex and noncommuting terms; uncached assembly provides the reference.

## Numerical / algorithmic content

Covers both cache insertion orders, changed coefficients, repeated hits, output arity, crystal full/Zeeman decomposition, and ordinary spin/mode cache controls.

## Syntax

`result=test_giant_cache()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
