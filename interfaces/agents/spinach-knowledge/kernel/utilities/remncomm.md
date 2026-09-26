# kernel/utilities/remncomm.m

- Signature: `A=remncomm(A,EvecB,EvalB)`

## Purpose

Removes from the Hermitian operator A the part that does not com- mute with the Hermitian operator B. Syntax: C=remncomm(A,EvecB,EvalB)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
- EvalB - a column vector containing the eigenvalues
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
