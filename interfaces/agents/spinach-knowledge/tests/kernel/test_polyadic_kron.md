# tests/kernel/test_polyadic_kron.m

- Signature: `result=test_polyadic_kron()`

## Purpose

Affixed polyadic tensor products and spin-space flow lifting.

## Physical / mathematical content

Kronecker products preserve every matrix factor; scalar and voxel descriptions of the same velocity must give the same spatial-spin generator.

## Numerical / algorithmic content

Complex rectangular and multiple-affix cases cover both operand orders, full/action references, uniform voxel flow, and divergence-free shear.

## Syntax

`result=test_polyadic_kron()`

## Parameters / inputs

None. The test constructs its own bounded physical fixtures.

## Outputs

`result` is the regression record of checks, messages, and failures; the test runner determines its final status.

## Header notes

The regression is registered in `test_manifest` and is not an optimisation or performance benchmark.
