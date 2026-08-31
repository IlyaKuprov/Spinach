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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Remove tiny negative round-off artefacts; implemented by `s=real(s(:))`.
- Lines 28-29: Check consistency; implemented by `grumble(s,tol)`.
- Lines 31-32: Project any remaining tiny negative round-off to zero; implemented by `s=max(s,0)`.
- Lines 34-35: Find the cutting point; implemented by `x=cumsum(s(end:-1:1).^2)`.
- Lines 38-39: Treat the zero case; implemented by `if isempty(k)`.

### Control flow inferred from the code

- Line 39: conditional branch on `isempty(k)`.

### Key state/data transformations

- Lines 25: computes `s` using `s=real(s(:))`.
- Lines 26: computes `s(abs(s)<1e-12)` using `s(abs(s)<1e-12)=0`.
- Lines 35: computes `x` using `x=cumsum(s(end:-1:1).^2)`.
- Lines 40: computes `r` using `r=0`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(s,tol)`.
  - Representative operation: `if (~isnumeric(tol))||(~isreal(tol))||(~isscalar(tol))||(tol<0)`.
  - Representative operation: `error('tol must be a non-negative real scalar.')`.

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
