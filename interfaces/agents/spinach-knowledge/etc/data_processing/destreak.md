# etc/data_processing/destreak.m

- Signature: `spectrum=destreak(spectrum)`

## Purpose

Reduces streak artefacts in 2D and 3D NMR spectra. Edges of the input spectrum must be free of genuine signals. Cell arrays and structures are processed recursively. Syntax: spectrum=destreak(spectrum)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- spectrum -a 2D or a 3D array, or a cell
- array, or a structure thereof

## Outputs

- spectrum -a 2D or a 3D array, or a cell
- array, or a structure thereof
- The function works by subtracting the kronecker propduct of edge
- lines from the spectrum matrix.

## Implementation structure

- Reduces streak artefacts in 2D and 3D NMR spectra. Edges of the
- input spectrum must be free of genuine signals. Cell arrays and
- structures are processed recursively. Syntax:
- spectrum=destreak(spectrum)
- spectrum -a 2D or a 3D array, or a cell
- array, or a structure thereof
- The function works by subtracting the kronecker propduct of edge
- lines from the spectrum matrix.
- Process structures and cell arrays recursively
- Get the field names
- Loop over structure elements
- Loop over field names
