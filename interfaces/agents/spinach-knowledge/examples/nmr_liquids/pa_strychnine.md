# examples/nmr_liquids/pa_strychnine.m

- Signature: `pa_strychnine()`

## Purpose

Simulates the ¹H NMR spectrum of strychnine, including a Redfield-superoperator model of line widths. The source estimates a calculation time of seconds.

## Physical / mathematical content

- Uses a Redfield relaxation superoperator to model line broadening in the liquid-state spin system.

## Numerical / algorithmic content

- Acquires the liquid-state FID, applies exponential apodisation, Fourier-transforms it, and plots the real spectrum.

## Implementation structure

- Configures the strychnine spin system and Redfield relaxation, then performs acquisition and signal processing.
