# kernel/utilities/kronm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/kronm.m`
- Signature: `x=kronm(Q,x)`
- Total lines: 111

## Purpose

Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening Kronecker products. Syntax: x=kronm(Q,x)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- Q -cell array of Kronecker terms
- x -a vector or a matrix of appropriate dimension
- Output:
- x -a vector or a matrix of appropriate dimension

## Implementation structure

- Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*x without opening
- Kronecker products. Syntax:
- x=kronm(Q,x)
- Q - cell array of Kronecker terms
- x - a vector or a matrix of appropriate dimension
- Output:
- Check consistency
- Number of matrices in Q
- Number of columns in x
- Row and column counts in Q
- Dimension map for x
- Reshape into the map

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `row_dims()`, `col_dims()`, `x_dims()`, `dims()`, `ipermute()`, `iscell()`, `ismatrix()`.
