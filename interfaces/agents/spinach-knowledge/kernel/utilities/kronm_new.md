# kernel/utilities/kronm_new.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/kronm_new.m`
- Signature: `M=kronm_new(Q,M)`
- Total lines: 71

## Purpose

Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening Kronecker products. Syntax: M=kronm(Q,M)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- Q -cell array of Kronecker terms
- M -a vector or a matrix of appropriate dimension
- Output:
- M -a vector or a matrix of appropriate dimension

## Implementation structure

- Calculates (Q{1}(x)Q{2}(x)...(x)Q{n})*M without opening
- Kronecker products. Syntax:
- M=kronm(Q,M)
- Q - cell array of Kronecker terms
- M - a vector or a matrix of appropriate dimension
- Output:
- Check consistency
- Dimension statistics
- Row and column counts in Q
- Fold up implicit dimensions of M
- Run the products
- Contract each implicit dimension

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `row_dims()`, `col_dims()`, `tensorprod()`, `iscell()`, `ismatrix()`.
