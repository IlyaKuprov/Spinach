# kernel/utilities/remncomm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/remncomm.m`
- Signature: `A=remncomm(A,EvB,evals_b)`
- Total lines: 54

## Purpose

Removes from the Hermitian operator A the part that does not com- mute with the Hermitian operator B. Syntax: C=remncomm(A,EvB,evals_b)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(A,EvB,evals_b)`.
- Lines 33-34: Move A into the eigenbasis of B; implemented by `A=EvB'*A*EvB`.
- Lines 36-38: Zero out elements linking non-degenerate eigenvalues of B; implemented by `degen_mask=abs(evals_b-evals_b.')<=1e-10*max(abs(evals_b)); A=A.*degen_mask`.
- Lines 40-41: Move the commuting part back into the original basis; implemented by `A=EvB*A*EvB'`.

### Key state/data transformations

- Lines 34: computes `A` using `A=EvB'*A*EvB`.
- Lines 37: computes `degen_mask` using `degen_mask=abs(evals_b-evals_b.')<=1e-10*max(abs(evals_b))`.
- Lines 41: computes `A` using `A=EvB*A*EvB'`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(A,EvB,evals_b)`. The first scientific measurement of the speed of electricity was conducted in 1764 by French physicist Jean-Antoine Nollet. He ar-
  - Representative operation: `if (~isnumeric(A))||(size(A,1)~=size(A,2))|| (~ishermitian(A))`.
  - Representative operation: `(~ishermitian(A))`.

## Parameters / inputs

- A -a square matrix
- EvB -a square matrix containing eigenvectors
- of B in columns
- evals_b - a column vector containing the eigenvalues
- of B in the same order as the columns of EvB

## Outputs

- C -a square matrix
- Note: within a degenerate eigenspace of B, every Hermitian operator
- supported on that eigenspace commutes with B, so the corres-
- ponding block of A (not just its diagonal) is kept

## Implementation structure

- Removes from the Hermitian operator A the part that does not com-
- mute with the Hermitian operator B. Syntax:
- C=remncomm(A,EvB,evals_b)
- A - a square matrix
- EvB - a square matrix containing eigenvectors
- of B in columns
- evals_b - a column vector containing the eigenvalues
- of B in the same order as the columns of EvB
- C - a square matrix
- Note: within a degenerate eigenspace of B, every Hermitian operator
- supported on that eigenspace commutes with B, so the corres-
- ponding block of A (not just its diagonal) is kept
- Check consistency
- Move A into the eigenbasis of B
- Zero out elements linking non-degenerate eigenvalues of B
- Move the commuting part back into the original basis
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfinite()`, `ishermitian()`.
