# kernel/utilities/frob_chop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/frob_chop.m`
- Signature: `r=frob_chop(s,tol)`
- Total lines: 61

## Purpose

Truncates SVD decomposition to the user-specified threshold in the Frobenius norm. Syntax: r=frob_chop(s,tol)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- s -a vector of singular values for a matrix,
- in descending order
- tol -truncation threshold

## Outputs

- r -the number of singular values to keep

## Implementation structure

- Truncates SVD decomposition to the user-specified threshold
- in the Frobenius norm. Syntax:
- r=frob_chop(s,tol)
- s -a vector of singular values for a matrix,
- in descending order
- tol -truncation threshold
- r -the number of singular values to keep
- Remove tiny negative round-off artefacts
- Check consistency
- Project any remaining tiny negative round-off to zero
- Find the cutting point
- Treat the zero case

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cumsum()`, `isscalar()`, `isvector()`, `any()`.
