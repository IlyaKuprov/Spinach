# kernel/pulses/pmlg5.m

- Signature: `phi=pmlg5(n)`

## Purpose

PMLG5 phase sequence as described in the paper by Vinogradova, Madhu and Vega (https://doi.org/10.1016/S0009-2614(99)01174-4).

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Syntax

```matlab
phi=spinal(n)
```

## Parameters / inputs

- n -a positive integer number

## Outputs

- phi -the phase of the n-th pulse in
- PMLG sequence, radians

## Implementation structure

- PMLG5 phase sequence as described in the paper by Vinogradova,
- Madhu and Vega (https://doi.org/10.1016/S0009-2614(99)01174-4).
- phi=spinal(n)
- n -a positive integer number
- phi -the phase of the n-th pulse in
- PMLG sequence, radians
- Check consistency
- PMLG5 phase sequence
- Loop correctly over
- Consistency enforcement
- One man's crappy software is another
- man's full time job.
