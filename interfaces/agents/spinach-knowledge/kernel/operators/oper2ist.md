# kernel/operators/oper2ist.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/oper2ist.m`
- Signature: `[states,coeffs]=oper2ist(A)`
- Total lines: 55

## Purpose

Irreducible spherical tensor operator expansion of a user- specified square matrix. Syntax: [states,coeffs]=oper2ist(A)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(A)`.
- Lines 30-31: Spherical tensors; implemented by `T=irr_sph_ten(size(A,1))`.
- Lines 33-34: Get expansion coefficients and states; implemented by `coeffs=cellfun(@(X)full(hdot(X,A)/hdot(X,X)),T)`.
- Lines 37-38: Drop negligible states; implemented by `idx=(abs(coeffs)>10*eps('double'))`.

### Key state/data transformations

- Lines 31: computes `T` using `T=irr_sph_ten(size(A,1))`.
- Lines 34: computes `coeffs` using `coeffs=cellfun(@(X)full(hdot(X,A)/hdot(X,X)),T)`.
- Lines 35: computes `states` using `states=transpose(0:(numel(T)-1))`.
- Lines 38: computes `idx` using `idx=(abs(coeffs)>10*eps('double'))`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(A)`. The obedient always think of themselves as virtuous rather than cowardly.
  - Representative operation: `if (~isnumeric(A))||(~ismatrix(A))|| (size(A,1)~=size(A,2))`.
  - Representative operation: `(size(A,1)~=size(A,2))`.

## Parameters / inputs

- A -a square matrix

## Outputs

- states -states, in the Spinach IST basis index-
- ing convention, that contribute to the
- operator in question; use lin2lm() fun-
- ction to convert to L,M spherical tens-
- or indices
- coeffs -coefficients with which the ISTs enter
- the linear combination

## Implementation structure

- Irreducible spherical tensor operator expansion of a user-
- specified square matrix. Syntax:
- [states,coeffs]=oper2ist(A)
- A -a square matrix
- states -states, in the Spinach IST basis index-
- ing convention, that contribute to the
- operator in question; use lin2lm() fun-
- ction to convert to L,M spherical tens-
- or indices
- coeffs -coefficients with which the ISTs enter
- the linear combination
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `irr_sph_ten()`, `cellfun()`, `hdot()`, `transpose()`, `eps()`, `coeffs()`, `states()`, `ismatrix()`.
