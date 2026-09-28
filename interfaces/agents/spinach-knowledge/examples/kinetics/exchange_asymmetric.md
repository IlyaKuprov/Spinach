# examples/kinetics/exchange_asymmetric.m

- Signature: `exchange_asymmetric()`

## Purpose

Simulates a two-site asymmetric chemical-exchange pattern for two `1H` spin environments at 14.1 T. Their scalar offsets are 0 and 3, and the exchange-rate matrix is `[-500 2000; 500 -2000]`; the specified concentration weights are `[2000 500]`. The source lists a calculation time of seconds.

## Physical / mathematical content

The two exchanging environments have unequal forward and reverse rates, so their state interconversion is asymmetric. The example generates a liquid-state NMR signal from this exchange system.

## Numerical / algorithmic content

Uses the sphten-liouv formalism with no basis approximation, acquires the signal, applies exponential apodisation with parameter 6, and Fourier-transforms the zero-filled FID.

## Implementation structure

Specifies the two-spin exchange system, constructs the Spinach basis, sets the acquisition parameters, simulates the FID, and plots its Fourier-transformed spectrum.
