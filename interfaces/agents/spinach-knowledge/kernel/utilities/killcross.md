# kernel/utilities/killcross.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/killcross.m`
- Signature: `M=killcross(M,f1idx,f2idx)`
- Total lines: 55

## Purpose

Zeroes the specified rows and columns of a matrix. Syntax: M=killcross(M,f1idx,f2idx)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- M -a matrix
- f1idx -numbers of the columns that
- should be zeroed
- f2idx -numbers of the rows that
- should be zeroed

## Outputs

- M -a matrix

## Implementation structure

- Zeroes the specified rows and columns of a matrix. Syntax:
- M=killcross(M,f1idx,f2idx)
- M -a matrix
- f1idx -numbers of the columns that
- should be zeroed
- f2idx -numbers of the rows that
- Check consistency
- Wipe the indices
- Consistency enforcement
- A narcissist is someone better-looking than you are.
- Gore Vidal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismatrix()`, `any()`, `isequal()`.
