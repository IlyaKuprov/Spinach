# interfaces/nmrpipe/fid2ascii.m

- Signature: `fid2ascii(file_name,fid)` (both inputs are required; there are no optional arguments).

## Purpose and data contract

Writes numeric Spinach FID data to a headerless, whitespace-separated ASCII table. The coordinate columns are one-based point indices, not physical times or frequencies; sweep widths and other acquisition metadata are not written and must be supplied to the downstream program.

- A vector is written as two columns: real samples at indices `1:N`, followed by imaginary samples at indices `N+1:2*N`. The imaginary block is written even when the input is real.
- A matrix is written as three columns: direct-dimension index, second-dimension index, and value. For each second-dimension index, all real values are written first, then imaginary values with the direct index shifted by the direct-dimension length.
- A three-dimensional array is written as four columns: first-, second-, and third-dimension indices, then value. The real block precedes the imaginary block for each pair of higher-dimensional indices, with the first index shifted by its dimension length for the imaginary block.

Values use the `%12.8E` numeric format. This is an ASCII point-index/value exchange format; it does not attach time, sweep-width, or unit metadata. Structure-valued FIDs are not accepted as a whole and must be exported field by field.

## Inputs and guardrails

- `file_name`: destination passed to `fopen(file_name,'w')`; the function does not validate its type or check the returned file identifier.
- `fid`: any numeric vector, matrix, or 3-D numeric array. The only explicit input guard is `isnumeric(fid)`; values are not checked for finiteness, and no shape/size consistency checks are made. Other dimensionalities raise `unsupported data dimensionality.`

The function writes the file and returns no value. The existing author contact is ilya.kuprov@weizmann.ac.il.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/nmrpipe/fid2ascii.m)
- [Spin Dynamics Wiki: fid2ascii.m](https://spindynamics.org/wiki/index.php?title=fid2ascii.m)
