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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(M,f1idx,f2idx)`.
- Lines 28-29: Wipe the indices; implemented by `M(f2idx,:)=0; M(:,f1idx)=0`.

### Key state/data transformations

- Lines 29: computes `M(f2idx,:)` using `M(f2idx,:)=0; M(:,f1idx)=0`.

### Local helper functions

- Line 34: `grumble()` — `function grumble(M,f1idx,f2idx)`.
  - Representative operation: `if (~isnumeric(M))||(~ismatrix(M))`.
  - Representative operation: `error('M must be a matrix.')`.

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
