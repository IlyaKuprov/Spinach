# tests/kernel/test_shaped_pulse_rotation.m

- Signature: `result=test_shaped_pulse_rotation()`

## Purpose

Verifies that a one-slice Cartesian X pulse reproduces the hard-pulse rotation limit.

## Physical / mathematical content

A rectangular X pulse with amplitude `1 rad/s` and duration `pi s` has flip angle `pi` radians, so it must invert `Lz` to `-Lz`.

## Numerical / algorithmic content

The test uses a one-proton Hilbert-space system with zero drift and controls `{Lx,Ly}`. It applies one slice with amplitudes `{1,0}` and duration `pi` using `shaped_pulse_xy` with the `expm-pwc` method, then compares the result with `-Lz` at absolute and relative tolerances of `1e-14`.

## Outputs

`result` is the regression-test record, including the rotation comparison and explanatory messages.

## Implementation structure

Constructs the spin system and Cartesian controls, applies the shaped pulse, and checks the expected state inversion.

## Header notes

The source header credits Ilya Kuprov (`ilya.kuprov@weizmann.ac.il`).
