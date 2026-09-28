# examples/nmr_liquids/pa_menthol.m

- Signature: `pa_menthol()`

## Purpose

Simulates the menthol NMR spectrum and the effects of poor Z1 and Z2 magnet shims, using spectrum information from Damien Jeannerat. The source estimates a calculation time of minutes.

## Physical / mathematical content

- Represents the menthol proton spin system with scalar couplings and models the specified Z1/Z2 shim errors.

## Numerical / algorithmic content

- Builds a scalar-coupling Liouville-space model and simulates a liquid-state FID.
- Applies Gaussian apodisation and the bad-shim effects, Fourier-transforms the signal, and plots the spectrum.

## Implementation structure

- Sets up the menthol spin system, acquisition, and shim parameters; computes and processes the FID before plotting.
