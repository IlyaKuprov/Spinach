# tests/kernel/test_srsk_thermal.m

- Signature: `result=test_srsk_thermal()`

## Purpose

Tests once-only thermalisation of additive SRSK relaxation and once-only addition of bosonic mode dissipation.

## Physical / mathematical content

A rapidly relaxing `14N` source is coupled to `1H` with positive, negative, or zero scalar coupling. An oriented quadrupolar interaction supplies a complex, noncommuting thermalisation Hamiltonian. A spectator cavity distinguishes unital dephasing from non-unital amplitude damping at finite temperature and checks trace conservation.

## Numerical / algorithmic content

Compares zero-destination SRSK rates with explicit rate augmentation, and IME and DiBari thermalisation with once-only references across coupling signs and supported retention policies. Spin-boson cases cover damping, dephasing, both, and neither, at nonzero and zero scalar coupling, against spin-only thermalisation plus one original-temperature mode dissipator and the no-SRSK path with explicitly augmented rates.

## Syntax

`result=test_srsk_thermal()`

## Parameters / inputs

None. The test constructs its own fixtures.

## Outputs

`result` contains regression checks for rates, retention, and equilibrium.