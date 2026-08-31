# kernel/conventions/transforms/stev2sph.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/stev2sph.m`
- Signature: `Bkq=stev2sph(k,Bkq)`
- Total lines: 79

## Purpose

Transforms the coefficients in front of Stevens operators, as produced by stevens.m, into the coefficients before the irredu- cible spherical tensor operators, as produced by irr_sph_ten.m function. Works up to 6th spherical rank. Source:

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(k,Bkq)`.
- Lines 37-38: Catalog the stupid scaling factors; implemented by `a{1}=[1/sqrt(2) 1 1/sqrt(2)]'`.
- Lines 45-46: Form the transformation matrix diagonal; implemented by `criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1) ].*a{k}`.
- Lines 48-49: Form the transformation matrix antidiagonal; implemented by `cross=[+1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`.
- Lines 51-52: Form the transformation matrix; implemented by `A=diag(criss)+fliplr(diag(cross))`.
- Lines 54-55: Transform the coefficients; implemented by `Bkq=transpose(Bkq'*A)`.

### Key state/data transformations

- Lines 38: computes `a{1}` using `a{1}=[1/sqrt(2) 1 1/sqrt(2)]'`.
- Lines 39: computes `a{2}` using `a{2}=[1 1/2 sqrt(6) 1/2 1]'`.
- Lines 40: computes `a{3}` using `a{3}=[sqrt(2) 1/sqrt(3) sqrt(10/3) sqrt(10) sqrt(10/3) 1/sqrt(3) sqrt(2)]'`.
- Lines 41: computes `a{4}` using `a{4}=[2 1/sqrt(2) sqrt(7) sqrt(7/2) 2*sqrt(70) sqrt(7/2) sqrt(7) 1/sqrt(2) 2]'`.
- Lines 42: computes `a{5}` using `a{5}=[2*sqrt(2) 2/sqrt(5) 6*sqrt(2/5) sqrt(3/5) 2*sqrt(21/5) 6*sqrt(14) 2*sqrt(21/5) sqrt(3/5) 6*sqrt(2/5) 2/sqrt(5) 2*sqrt(2)]'`.
- Lines 43: computes `a{6}` using `a{6}=[4 2/sqrt(3) 4*sqrt(11/6) 2*sqrt(11/5) 4*sqrt(11/5) sqrt(22) 4*sqrt(231) sqrt(22) 4*sqrt(11/5) 2*sqrt(11/5) 4*sqrt(11/6) 2/sqrt(3) 4]'`.
- Lines 46: computes `criss` using `criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1) ].*a{k}`.
- Lines 49: computes `cross` using `cross=[+1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`.
- Lines 52: computes `A` using `A=diag(criss)+fliplr(diag(cross))`.
- Lines 55: computes `Bkq` using `Bkq=transpose(Bkq'*A)`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(k,Bkq)`. K.W.H. Stevens has done a great disservice to Magnetic Resonance by
  - Representative operation: `if (~isnumeric(k))||(~isreal(k))||(~isfinite(k))|| (~isscalar(k))||(mod(k,1)~=0)||(k<1)||(k>6)`.
  - Representative operation: `(~isscalar(k))||(mod(k,1)~=0)||(k<1)||(k>6)`.

## Syntax

```matlab
Bkq=stev2sph(k,Bkq)
```

## Parameters / inputs

- k -the spherical rank in question
- Bkq -a column of 2k+1 real coefficients
- in front of Stevens operators, in
- increasing order of projections

## Outputs

- Bkq -a column of 2k+1 complex coefficients
- in front of irreducible spherical
- tensor operators, in decreasing order
- of projections

## Implementation structure

- Transforms the coefficients in front of Stevens operators, as
- produced by stevens.m, into the coefficients before the irredu-
- cible spherical tensor operators, as produced by irr_sph_ten.m
- function. Works up to 6th spherical rank. Source:
- Bkq=stev2sph(k,Bkq)
- k -the spherical rank in question
- Bkq -a column of 2k+1 real coefficients
- in front of Stevens operators, in
- increasing order of projections
- Bkq -a column of 2k+1 complex coefficients
- in front of irreducible spherical
- tensor operators, in decreasing order

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fliplr()`, `transpose()`, `isscalar()`, `any()`, `iscolumn()`.
