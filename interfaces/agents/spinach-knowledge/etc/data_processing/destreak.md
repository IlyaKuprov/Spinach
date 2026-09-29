# etc/data_processing/destreak.m

- MATLAB implementation: [etc/data_processing/destreak.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/data_processing/destreak.m)

- Signature: `spectrum=destreak(spectrum)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=Destreak.m)
- Source author: ilya.kuprov@weizmann.ac.il

## Purpose

Reduce streak artefacts in 2D and 3D NMR spectra by subtracting repeated edge contributions. The first-column/first-plane edge and the other first-index edges used by the correction must be free of genuine signal; otherwise genuine data there will also be subtracted.

## Call and input

`spectrum=destreak(spectrum)`

- A numeric 2D or 3D array is processed directly.
- A cell array or structure is traversed recursively; each contained value is passed back to `destreak`, so its eventual numeric leaves must be supported arrays.
- Numeric vectors are rejected as one-dimensional spectra, and numeric arrays with dimensionality other than 2 or 3 produce an unsupported-dimensionality error.
- There are no scale or unit parameters; array values retain their units and the output has the same container and array shape.

## Mechanism

For a 2D matrix, the routine subtracts the first column replicated across columns, then subtracts the first row replicated across rows. For 3D data it successively subtracts the first-index face along each of the three dimensions, replicated over that dimension. This is a direct edge-based correction, not an estimate of which edge points are artefacts: the signal-free-edge condition is essential.

## Output and limitations

The returned `spectrum` contains the corrected numeric data in the original array, cell, or struct layout. Only 2D and 3D numeric leaves are supported; genuine edge signal is not protected from subtraction. The cell-recursion loop is written as `parfor`; the source does not specify a separate runtime or speed guarantee.
