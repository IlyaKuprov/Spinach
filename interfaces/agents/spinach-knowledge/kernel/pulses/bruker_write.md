# kernel/pulses/bruker_write.m

- Signature: `bruker_write(X,Y,dt,file_name)`

## Purpose

Writes a shaped pulse as a Bruker JCAMP text file for use in TopSpin. The routine converts the Cartesian pulse components to amplitude and phase, wraps phase into one turn and expresses it in degrees, and scales the amplitudes to Bruker's 0–100 range when the maximum amplitude is positive.

## Inputs

- `X` — real numeric column vector of pulse components in Hz.
- `Y` — real numeric column vector of pulse components in Hz, with the same number of elements as `X`.
- `dt` — positive real scalar duration of each time slice, in seconds.
- `file_name` — output filename as a character vector.

## Output and file contents

The function writes an ASCII Bruker shape file. Its header includes the pulse duration in microseconds, the number of points, amplitude and phase ranges, and the JCAMP shape metadata. The data section contains one amplitude/phase pair per pulse point, followed by `##END`.

## Reference

[Spin Dynamics Wiki: `bruker_write.m`](https://spindynamics.org/wiki/index.php?title=bruker_write.m)
