# kernel/pulses/spinal.m

MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/spinal.m
Source Wiki page: https://spindynamics.org/wiki/index.php?title=spinal.m

- Signature: `phi=spinal(n)`

## Purpose

Returns the phase, in radians, of pulse `n` in the SPINAL sequence described by Fung, Khitrin, and Ermolaev.

## Sequence and indexing

The implementation stores 64 phase entries in degrees, selects entry `mod(n-1,64)+1`, then converts that entry to radians. Indexing starts at 1, and the stored phase pattern repeats every 64 pulses; for example, the first entry is 10 degrees. The function returns only a phase value, not an RF amplitude or pulse duration.

## Parameters / inputs

- `n` - positive integer scalar identifying the pulse in the sequence

## Output

- `phi` - phase of the `n`th pulse in the SPINAL sequence, radians

## Reference

Fung, Khitrin, and Ermolaev, https://doi.org/10.1006/jmre.1999.1896
