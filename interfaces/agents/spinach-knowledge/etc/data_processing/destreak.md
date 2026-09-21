# etc/data_processing/destreak.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/data_processing/destreak.m`
- Signature: `spectrum=destreak(spectrum)`
- Total lines: 105

## Purpose

Reduces streak artefacts in 2D and 3D NMR spectra. Edges of the input spectrum must be free of genuine signals. Cell arrays and structures are processed recursively. Syntax: spectrum=destreak(spectrum)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isstruct()`, `fieldnames()`, `spectrum()`, `iscell()`, `grumble()`, `ndims()`, `isvector()`.
