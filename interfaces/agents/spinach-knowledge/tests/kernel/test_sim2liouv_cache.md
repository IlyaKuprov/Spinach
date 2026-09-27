# tests/kernel/test_sim2liouv_cache.m

- Signature: `result=test_sim2liouv_cache()`

## Purpose

Representation-specific caches after Hilbert-to-Liouville conversion.

## Physical / mathematical content

An operator or Hamiltonian must retain the dimensions and action of its requested representation.

## Numerical / algorithmic content

Checks both insertion orders, separate and combined caches, non-Hermitian operators, no-op formalisms, absent/disabled metadata, and cached soft-pulse acquisition.

## Syntax

`result=test_sim2liouv_cache()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
