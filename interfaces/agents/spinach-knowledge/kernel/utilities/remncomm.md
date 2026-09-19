# kernel/utilities/remncomm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/remncomm.m`
- Signature: `A=remncomm(A,EvecB,EvalB)`
- Total lines: 73

## Purpose

Removes from the Hermitian operator A the part that does not com- mute with the Hermitian operator B. Syntax: C=remncomm(A,EvecB,EvalB)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(A,EvecB,EvalB)`.
- Lines 33-34: Move A into the eigenbasis of B; implemented by `A=EvecB'*A*EvecB`.
- Lines 36-38: Zero out elements linking eigenvalues of B that differ by more than eigensolver roundoff; implemented by `degen_mask=abs(EvalB-EvalB.')<=100*numel(EvalB)*eps(max(EvalB)-min(EvalB)); A=A.*degen_mask`.
- Lines 40-41: Move the commuting part back into the original basis; implemented by `A=EvecB*A*EvecB'`.

### Key state/data transformations

- Lines 34: computes `A` using `A=EvecB'*A*EvecB`.
- Lines 37: computes `degen_mask` using `degen_mask=abs(EvalB-EvalB.')<=100*numel(EvalB)*eps(max(EvalB)-min(EvalB))`.
- Lines 41: computes `A` using `A=EvecB*A*EvecB'`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(A,EvecB,EvalB)`. The first scientific measurement of the speed of electricity was conducted in 1764 by French physicist Jean-Antoine Nollet. He ar-
  - Representative operation: `if (~isnumeric(A))||(size(A,1)~=size(A,2))|| (~ishermitian(A))`.
  - Representative operation: `(~ishermitian(A))`.

## Parameters / inputs

- A -a square matrix
- EvecB -a square matrix containing eigenvectors
- of B in columns
- EvalB -a column vector containing the eigenvalues
- of B in the same order as the columns of EvecB

## Outputs

- C -a square matrix
- Note: within a degenerate eigenspace of B, every Hermitian operator
- supported on that eigenspace commutes with B, so the corres-
- ponding block of A (not just its diagonal) is kept

## Implementation structure

- Removes from the Hermitian operator A the part that does not com-
- mute with the Hermitian operator B. Syntax:
- C=remncomm(A,EvecB,EvalB)
- A - a square matrix
- EvecB - a square matrix containing eigenvectors
- of B in columns
- EvalB -a column vector containing the eigenvalues
- of B in the same order as the columns of EvecB
- C - a square matrix
- Note: within a degenerate eigenspace of B, every Hermitian operator
- supported on that eigenspace commutes with B, so the corres-
- ponding block of A (not just its diagonal) is kept
- Check consistency
- Move A into the eigenbasis of B
- Zero out elements linking eigenvalues of B that differ by more than eigensolver roundoff
- Move the commuting part back into the original basis
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfinite()`, `ishermitian()`.
