# kernel/utilities/killdiag.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/killdiag.m`
- Signature: `spec=killdiag(spec,brush_dim)`
- Total lines: 63

## Purpose

Zeroes out the diagonal of a 2D spectrum using the brush with the specified dimensions. Syntax: spec=killdiag(spec,brush_dim)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(spec,brush_dim)`.
- Lines 26-27: Loop over the column index; implemented by `for n=1:size(spec,2)`.
- Lines 29-30: Find the row index; implemented by `k=n*size(spec,1)/size(spec,2)`.
- Lines 32-33: Find the row index extents; implemented by `k=floor(k-brush_dim/2):ceil(k+brush_dim/2)`.
- Lines 35-36: Avoid array boundaries; implemented by `k(k<1)=[]; k(k>size(spec,1))=[]`.
- Lines 38-39: Zero the elements; implemented by `spec(k,n)=0`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=1:size(spec,2)`.

### Key state/data transformations

- Lines 30: computes `k` using `k=n*size(spec,1)/size(spec,2)`.
- Lines 36: computes `k(k<1)` using `k(k<1)=[]; k(k>size(spec,1))=[]`.
- Lines 39: computes `spec(k,n)` using `spec(k,n)=0`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(spec,brush_dim)`.
  - Representative operation: `if (~isnumeric(spec))||(~ismatrix(spec))`.
  - Representative operation: `error('spec must be a matrix.')`.

## Parameters / inputs

- spec -2D matrix representing a spectrum
- brush_dim -the width of the band to zero out
- around the diagonal, points

## Outputs

- spec -2D matrix representing a spectrum

## Implementation structure

- Zeroes out the diagonal of a 2D spectrum using the brush
- with the specified dimensions. Syntax:
- spec=killdiag(spec,brush_dim)
- spec -2D matrix representing a spectrum
- brush_dim -the width of the band to zero out
- around the diagonal, points
- Check consistency
- Loop over the column index
- Find the row index
- Find the row index extents
- Avoid array boundaries
- Zero the elements

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spec()`, `ismatrix()`, `isscalar()`, `any()`.
