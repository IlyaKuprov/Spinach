# examples/kinetics/exchange_symmetric.m

- Signature: `exchange_symmetric()`

## Purpose

Simulates a two-site symmetric chemical-exchange pattern for two `1H` spin environments at 14.1 T. Their scalar offsets are 0 and 3, both exchange directions have rate 2000, and the concentration weights are `[1 1]`. The source lists a calculation time of seconds.

## Physical / mathematical content

The two environments interconvert at equal rates and have equal specified concentration weights. The example generates a liquid-state NMR signal from this exchange system.

## Numerical / algorithmic content

Uses the sphten-liouv formalism with no basis approximation, acquires the signal, applies exponential apodisation with parameter 6, and Fourier-transforms the zero-filled FID.

## Implementation structure

Specifies the two-spin exchange system, constructs the Spinach basis, sets the acquisition parameters, simulates the FID, and plots its Fourier-transformed spectrum.
