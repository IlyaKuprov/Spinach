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
