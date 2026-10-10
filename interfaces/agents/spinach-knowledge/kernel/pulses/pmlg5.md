# kernel/pulses/pmlg5.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/pmlg5.m
Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=pmlg5.m

## Purpose

Return the phase of a pulse in the PMLG5 phase cycle. This function returns one phase value, not a sampled RF waveform.

## Syntax

~~~matlab
phi=pmlg5(n)
~~~

## Input and output

- n: positive integer pulse index.
- phi: phase in radians.

## Implementation

The source stores a 20-value phase cycle in degrees, selects entry mod(n-1,20)+1, then multiplies it by pi/180. Thus indices beyond 20 repeat the cycle; the first entry is 339.22 degrees (converted to radians). The cited phase sequence is from Vinogradova, Madhu and Vega: https://doi.org/10.1016/S0009-2614(99)01174-4.
