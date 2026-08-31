# kernel/operators/oper2bm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/oper2bm.m`
- Signature: `[states,coeffs]=oper2bm(A)`
- Total lines: 71

## Purpose

Bosonic monomial operator expansion of a user-specified square matrix. Syntax: [states,coeffs]=oper2bm(A)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(A)`.
- Lines 30-31: Bosonic monomials; implemented by `BM=boson_mono(size(A,1))`.
- Lines 34-35: Compute the overlap matrix; implemented by `S=zeros(nlevels^2,nlevels^2)`.
- Lines 42-43: Get the expansion coefficients; implemented by `coeffs=zeros(nlevels^2,1)`.
- Lines 49-50: List the states; implemented by `states=transpose(0:(numel(BM)-1))`.
- Lines 52-53: Drop negligible states; implemented by `idx=(abs(coeffs)>10*eps('double'))`.

### Control flow inferred from the code

- Line 36: `for` loop over `n=1:nlevels^2`.
- Line 37: `for` loop over `k=1:nlevels^2`.
- Line 44: `for` loop over `n=1:nlevels^2`.

### Key state/data transformations

- Lines 31: computes `BM` using `BM=boson_mono(size(A,1))`.
- Lines 32: computes `nlevels` using `nlevels=size(A,1)`.
- Lines 35: computes `S` using `S=zeros(nlevels^2,nlevels^2)`.
- Lines 38: computes `S(n,k)` using `S(n,k)=hdot(BM{n},BM{k})`.
- Lines 43: computes `coeffs` using `coeffs=zeros(nlevels^2,1)`.
- Lines 45: computes `coeffs(n)` using `coeffs(n)=hdot(BM{n},A)`.
- Lines 50: computes `states` using `states=transpose(0:(numel(BM)-1))`.
- Lines 53: computes `idx` using `idx=(abs(coeffs)>10*eps('double'))`.

### Local helper functions

- Line 59: `grumble()` — `function grumble(A)`. The greatest trainers can teach important skills, but the best
  - Representative operation: `if (~isnumeric(A))||(~ismatrix(A))|| (size(A,1)~=size(A,2))`.
  - Representative operation: `(size(A,1)~=size(A,2))`.

## Parameters / inputs

- A -a square matrix

## Outputs

- states -states, in the Spinach BM basis index-
- ing convention, that contribute to the
- operator in question; use lin2kq() fun-
- ction to convert to K,Q bosonic monomi-
- al indices
- coeffs -coefficients with which the BMs enter
- the linear combination

## Implementation structure

- Bosonic monomial operator expansion of a user-specified
- square matrix. Syntax:
- [states,coeffs]=oper2bm(A)
- A -a square matrix
- states -states, in the Spinach BM basis index-
- ing convention, that contribute to the
- operator in question; use lin2kq() fun-
- ction to convert to K,Q bosonic monomi-
- al indices
- coeffs -coefficients with which the BMs enter
- the linear combination
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `boson_mono()`, `hdot()`, `coeffs()`, `transpose()`, `eps()`, `states()`, `ismatrix()`.
