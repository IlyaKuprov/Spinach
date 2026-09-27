# examples/nmr_liquids/pa_sucrose.m

- Signature: `pa_sucrose()`

## Purpose

Simulates the ¹H NMR spectrum of sucrose using magnetic parameters read from a DFT calculation and a Redfield relaxation superoperator. The source estimates a calculation time of seconds.

## Physical / mathematical content

- Uses DFT-derived magnetic parameters for the sucrose spin system and includes Redfield relaxation.

## Numerical / algorithmic content

- Builds the relaxation-enabled model, simulates the liquid-state FID, applies exponential apodisation, Fourier-transforms it, and plots the spectrum.

## Implementation structure

- Imports the magnetic parameters, configures the spin system and Redfield relaxation, then runs acquisition and signal processing.
