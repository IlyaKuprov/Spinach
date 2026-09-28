# examples/fundamentals/strong_coupling.m

- Signature: `strong_coupling()`

## Purpose

A garden variety strongly coupled two-spin system.

## Physical / mathematical content

Models two 1H spins at 5.9 T with scalar shifts 1.0 and 1.5 and a scalar coupling of 7.0, retaining the full spin dynamics in the spherical-tensor Liouville formalism.

## Numerical / algorithmic content

Simulates a liquid-state acquisition with 1024 points, a 300 Hz sweep and offset, applies exponential apodisation with parameter 10, zero-fills to 4096 points, and Fourier-transforms the FID for plotting.

## Implementation structure

Creates the system and basis, configures a 1H acquisition with raising-operator initial state and receiver, runs the liquid simulation, processes the FID, and plots the real spectrum.
