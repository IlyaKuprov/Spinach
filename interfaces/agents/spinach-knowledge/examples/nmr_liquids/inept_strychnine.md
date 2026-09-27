# examples/nmr_liquids/inept_strychnine.m

- Signature: `inept_strychnine()`

## Purpose

INEPT experiment on strychnine. Calculation time: minutes.

## Implementation

The strychnine spin system is simulated at 5.9 T and 298 K using an IK-2 scalar-coupling basis. The example generates 13C isotopomers and runs the INEPT sequence with a 150 Hz transfer coupling for the 13C/1H spin channels. It exponentially apodises and Fourier transforms each FID before plotting the carbon spectrum.
