# examples/kinetics/flux_asymmetric.m

- Signature: `flux_asymmetric()`

## Purpose

Simulates a two-site asymmetric intermolecular magnetization-flux problem for two `1H` environments at 14.1 T, with scalar offsets 0 and 3. The flux rates are 500 from site 1 to 2 and 2000 from site 2 to 1. The source lists a calculation time of seconds.

## Physical / mathematical content

The model sets unequal directional intermolecular flux rates and computes the resulting liquid-state NMR signal.

## Numerical / algorithmic content

Uses the sphten-liouv formalism with no basis approximation, acquires the signal, applies exponential apodisation with parameter 6, and Fourier-transforms the zero-filled FID.

## Implementation structure

Specifies the two-spin flux system, constructs the Spinach basis, sets acquisition parameters, simulates the FID, and plots its Fourier-transformed spectrum.
