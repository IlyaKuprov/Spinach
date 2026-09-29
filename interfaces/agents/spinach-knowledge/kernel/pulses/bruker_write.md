# kernel/pulses/bruker_write.m

[Source: `kernel/pulses/bruker_write.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/bruker_write.m)

- Signature: `bruker_write(X,Y,dt,file_name)`

## Purpose

Exports Cartesian RF samples as a Bruker JCAMP shape file for TopSpin. For each paired X/Y Cartesian sample, `cartesian2polar` gives amplitude and phase; phase is wrapped to one turn and written in degrees. If the maximum amplitude is positive, amplitudes are scaled so that maximum is 100. This is a relative shape scale, not an exported Hz amplitude. All-zero input stays zero.

## Inputs and sampling

- `X`, `Y` — equal-length real numeric column vectors of Cartesian RF components, in Hz.
- `dt` — positive real slice duration in seconds. With `N = numel(X)`, the header pulse length is `N*dt`, converted to microseconds; there is one amplitude/phase pair for each of the N samples.
- `file_name` — character vector naming the output text file; the documented convention is a `.txt` extension (the implementation checks character type, not the extension).

The function checks that X and Y are numeric real column vectors of equal length, dt is a positive real scalar, and the filename is a character vector. It does not emit a MATLAB output value.

## File output

The header identifies Bruker JCAMP-DX shape data and Spinach, records the current date and time, amplitude and phase extrema, pulse length, and point count, then declares `##XYPOINTS= (XY..XY..)`. Each following row contains one normalised amplitude and its phase in degrees, separated by a space; `##END` closes the file. The initial `writelines(lines,file_name)` writes the header to the destination (replacing an existing file under MATLAB's default write behaviour); the numeric pairs and terminator are appended. This operation therefore has a file-system side effect and overwrites an existing destination rather than adding a new pulse to it.

## Reference

[Spin Dynamics Wiki: `bruker_write.m`](https://spindynamics.org/wiki/index.php?title=bruker_write.m)
