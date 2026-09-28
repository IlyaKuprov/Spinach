# interfaces/nmrpipe/fid2ascii.m

- Signature: `fid2ascii(file_name,fid)`

## Purpose

Writes a Spinach-generated one-, two-, or three-dimensional free induction decay (FID) to an ASCII file. Structure arrays must be exported field by field.

## Inputs and format

- `file_name`: output filename (a character string).
- `fid`: numeric FID array.
- Each row contains point-number coordinates followed by a value. For each combination of higher-dimensional coordinates, real values are written first, then imaginary values; first-dimension point numbers for the imaginary values are offset by the length of that dimension.
- Coordinates are point numbers, not physical frequencies. Sweep width is not stored and must be supplied to the downstream program. Unsupported dimensionality raises an error.

## Output

The function writes the ASCII file and returns no value.

## Source

Contact: ilya.kuprov@weizmann.ac.il. [fid2ascii.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=fid2ascii.m)
