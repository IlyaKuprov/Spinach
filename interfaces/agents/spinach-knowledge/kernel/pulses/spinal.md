# kernel/pulses/spinal.m

- Signature: `phi=spinal(n)`

## Purpose

SPINAL phase sequences as described in the paper by Fung, Khitrin and Ermolaev (https://doi.org/10.1006/jmre.1999.1896). Syntax: phi=spinal(n)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- n -a positive integer number

## Outputs

- phi -the phase of the n-th pulse in
- SPINAL sequence, radians

## Implementation structure

- SPINAL phase sequences as described in the paper by Fung, Khitrin
- and Ermolaev (https://doi.org/10.1006/jmre.1999.1896). Syntax:
- phi=spinal(n)
- n -a positive integer number
- phi -the phase of the n-th pulse in
- SPINAL sequence, radians
- Check consistency
- Spinal phase sequence
- Loop correctly over
- Consistency enforcement
- The key to performance is elegance, not
- battalions of special cases.
