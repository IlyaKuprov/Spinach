# examples/nmr_liquids/inv_rec_strychnine.m

- Signature: `inv_rec_strychnine()`

## Purpose

1H inversion-recovery experiment on strychnine at 250 MHz. Calculation time: minutes.

## Implementation

The code builds the 1H strychnine system at 5.9 T with Redfield relaxation, the Di Bari equilibrium state, and an IK-2 scalar-coupling basis. It simulates inversion recovery over ten delays up to 1 s, applies exponential apodisation to the FIDs, and Fourier transforms them for plotting.
