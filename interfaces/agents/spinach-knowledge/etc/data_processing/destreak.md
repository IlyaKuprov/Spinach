# etc/data_processing/destreak.m

- Signature: `spectrum=destreak(spectrum)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=Destreak.m)

## Purpose

Reduce streak artefacts in multidimensional NMR spectra by subtracting signal-free edge contributions. The input edges used for this correction must not contain genuine signals.

## Method

For a numeric 2D array, the routine subtracts contributions repeated from the first column and first row. For a 3D array, it subtracts the corresponding first-plane/edge contributions along each dimension. Structs and cell arrays are traversed recursively; each contained value is passed back to `destreak`.

## Input and output

- `spectrum` — numeric 2D or 3D array, cell array, or struct containing values processed recursively. Numeric vectors are rejected.
- `spectrum` — corrected data, preserving the input container and array shape.
