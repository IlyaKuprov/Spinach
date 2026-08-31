# kernel/utilities/keep_rank.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/keep_rank.m`
- Signature: `A=keep_rank(A,nsvk)`
- Total lines: 50

## Purpose

Truncates the singular value decomposition at the specified rank and reassembles the matrix. Syntax: A=keep_rank(A,rank)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(A,nsvk)`.
- Lines 27-28: Run singular value decomposition; implemented by `[U,S,V]=svd(full(A))`.
- Lines 30-31: Truncate to the specified rank and rebuild; implemented by `A=U(:,1:nsvk)*S(1:nsvk,1:nsvk)*V(:,1:nsvk)'`.

### Key state/data transformations

- Lines 28: computes `[U,S,V]` using `[U,S,V]=svd(full(A))`.
- Lines 31: computes `A` using `A=U(:,1:nsvk)*S(1:nsvk,1:nsvk)*V(:,1:nsvk)'`.

### Local helper functions

- Line 36: `grumble()` — `function grumble(A,nsvk)`. Endure. In enduring, grow strong.
  - Representative operation: `if (~isnumeric(A))||(size(A,1)<=1)||(size(A,2)<=1)`.
  - Representative operation: `error('A must be a matrix.')`.

## Parameters / inputs

- A -real or complex matrix, will be
- converted to full if a sparse
- matrix is received
- nsvk -number of singular values to keep

## Outputs

- A -filtered matrix, returned as full

## Implementation structure

- Truncates the singular value decomposition at the specified rank
- and reassembles the matrix. Syntax:
- A=keep_rank(A,rank)
- A - real or complex matrix, will be
- converted to full if a sparse
- matrix is received
- nsvk - number of singular values to keep
- A - filtered matrix, returned as full
- Check consistency
- Run singular value decomposition
- Truncate to the specified rank and rebuild
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `any()`, `dim()`.
