# kernel/utilities/iseye.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/iseye.m`
- Signature: `verdict=iseye(M)`
- Total lines: 68

## Purpose

Returns true for unit matrices. The test is designed to be computationally affordable. Syntax: verdict=iseye(M)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- M -a matrix

## Outputs

- verdict -true or false

## Implementation structure

- Returns true for unit matrices. The test is designed to be
- computationally affordable. Syntax:
- verdict=iseye(M)
- M -a matrix
- verdict -true or false
- Check consistency
- Run the checks
- Not even square
- Not even diagonal
- Test vector
- Compare with unit
- Test failed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `false()`, `isdiag()`, `nnz()`, `true()`.
