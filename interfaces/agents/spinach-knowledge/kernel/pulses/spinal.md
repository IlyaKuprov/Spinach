# kernel/pulses/spinal.m

- Signature: `phi=spinal(n)`

## Purpose

Returns the phase, in radians, of pulse `n` in the SPINAL sequence described by Fung, Khitrin, and Ermolaev.

## Numerical / algorithmic content

The implementation stores a 64-entry phase sequence in degrees. It selects entry `mod(n-1,64)+1` and converts that phase to radians, so the pattern repeats every 64 pulses.

## Parameters / inputs

- `n` - positive integer scalar identifying the pulse in the sequence

## Outputs

- `phi` - phase of the `n`th pulse in the SPINAL sequence, radians

## Reference

Fung, Khitrin, and Ermolaev, https://doi.org/10.1006/jmre.1999.1896

Source Wiki page: https://spindynamics.org/wiki/index.php?title=spinal.m
