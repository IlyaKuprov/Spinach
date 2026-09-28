# tests/kernel/test_chemical_exchange_conservation.m

- Signature: `result=test_chemical_exchange_conservation()`

## Purpose

Tests conservation of total spin population in a closed, symmetric two-site chemical-exchange model.

## Physical / mathematical content

For a closed Markov exchange generator, conservation is checked by verifying that every column of `K` sums to zero.

## Numerical / algorithmic content

The test computes `K=kinetics(spin_system)` and compares its column sums with zero using absolute and relative tolerances of `1e-15`.

## Parameters / inputs

The model has two `1H` sites, `sys.magnet=14.1`, exchange rates `[-3 3; 3 -3]`, concentrations `[1 1]`, `sphten-liouv` formalism, and no basis approximation.

## Outputs

`result` contains the regression test result with explanatory messages.
