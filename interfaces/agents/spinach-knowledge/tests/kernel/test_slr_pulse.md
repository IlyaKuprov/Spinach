# tests/kernel/test_slr_pulse.m

- Signature: `result=test_slr_pulse()`

## Purpose

Tests Shinnar–Le Roux (SLR) selective-excitation pulse design, waveform consistency, excitation profile, production-path propagation, and input validation.

## Physical / mathematical content

The representative design uses 64 samples over `4 ms`, time-bandwidth product `4`, a `pi/2` flip angle, and passband and stopband ripple targets of `0.01`. Independent spin-half propagation checks that the `pi/2` waveform maps `Lz` to `-Ly`, and that a `pi/6` waveform maps `Lz` to `cos(pi/6)*Lz-sin(pi/6)*Ly`.

## Numerical / algorithmic content

The test checks waveform dimensions and finiteness, duration sum, and consistency between Cartesian and polar controls. A `16,385`-point Cayley–Klein frequency sweep checks unitarity error below `2e-12` and the target passband/stopband bounds. It also propagates the generated waveform through `shaped_pulse_xy` using `expm-pwc` and compares with `-Ly` at relative and absolute tolerances of `1e-10` each. Invalid designs are checked for rejection, including 63 samples, zero duration, ripple `1`, flip angle `pi`, and an infeasible 8-sample design with time-bandwidth product `0.1`.

## Outputs

`result` is the regression-test record with explanatory messages.

## Implementation structure

Generates an SLR waveform, checks its controls, propagates independent spin-half reference cases, evaluates the excitation profile, tests the production shaped-pulse path, and exercises invalid inputs.

## Header notes

The source header credits Ilya Kuprov (`ilya.kuprov@weizmann.ac.il`).
