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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Process structures and cell arrays recursively; implemented by `if isstruct(spectrum)`.
- Lines 29-30: Get the field names; implemented by `struct_fieldnames=fieldnames(spectrum)`.
- Lines 32-33: Loop over structure elements; implemented by `for m=1:numel(spectrum)`.
- Lines 35-36: Loop over field names; implemented by `for n=1:length(struct_fieldnames)`.
- Lines 38-39: Call itself for each field of each element; implemented by `spectrum(m).(struct_fieldnames{n})=destreak(spectrum(m).(struct_fieldnames{n}))`.
- Lines 45-46: Done; implemented by `return`.
- Lines 50-51: Loop over cells; implemented by `parfor n=1:numel(spectrum)`.
- Lines 53-54: Call itself for each cell; implemented by `spectrum{n}=destreak(spectrum{n})`.
- Lines 63-64: Check consistency; implemented by `grumble(spectrum)`.
- Lines 66-67: Decide problem dimensionality; implemented by `switch ndims(spectrum)`.
- Lines 71-72: Destreak the spectrum; implemented by `spectrum=spectrum-kron(spectrum(:,1),ones(1,size(spectrum,2)))`.
- Lines 77-78: Destreak the spectrum; implemented by `spectrum=spectrum-repmat(spectrum(:,1,:),1,size(spectrum,2),1)`.
- Lines 84-85: Complain and bomb out; implemented by `error('unsupported spectrum dimensionality.')`.

### Control flow inferred from the code

- Line 27: conditional branch on `isstruct(spectrum)`.
- Line 33: `for` loop over `m=1:numel(spectrum)`.
- Line 36: `for` loop over `n=1:length(struct_fieldnames)`.
- Line 51: `parfor` loop over `n=1:numel(spectrum)`.
- Line 67: dispatches on `ndims(spectrum)`; cases `2`, `3`.

### Key state/data transformations

- Lines 30: computes `struct_fieldnames` using `struct_fieldnames=fieldnames(spectrum)`.
- Lines 39: computes `spectrum(m).(struct_fieldnames{n})` using `spectrum(m).(struct_fieldnames{n})=destreak(spectrum(m).(struct_fieldnames{n}))`.
- Lines 54: computes `spectrum{n}` using `spectrum{n}=destreak(spectrum{n})`.
- Lines 72: computes `spectrum` using `spectrum=spectrum-kron(spectrum(:,1),ones(1,size(spectrum,2)))`.

### Local helper functions

- Line 92: `grumble()` — `function grumble(spectrum)`. If it flies, floats or fucks, you are better off renting it.
  - Representative operation: `if ~isnumeric(spectrum)`.
  - Representative operation: `error(['spectrum must be a numeric array, a cell array ' 'of numeric arrays or a structure with numeric arrays inside.'])`.

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
