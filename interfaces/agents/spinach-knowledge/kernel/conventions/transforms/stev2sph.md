# kernel/conventions/transforms/stev2sph.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/stev2sph.m`
- Signature: `Bkq=stev2sph(k,Bkq)`
- Total lines: 92

## Purpose

Transforms the coefficients in front of Stevens operators, as produced by stevens.m, into the coefficients before the irredu- cible spherical tensor operators, as produced by irr_sph_ten.m function. Works up to 12th spherical rank. Source for ranks up to 6:

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(k,Bkq)`.
- Lines 44-45: Catalog the stupid scaling factors; implemented by `a{1}=[1/sqrt(2) 1 1/sqrt(2)]'`.
- Lines 58-59: Form the transformation matrix diagonal; implemented by `criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1) ].*a{k}`.
- Lines 61-62: Form the transformation matrix antidiagonal; implemented by `cross=[+1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`.
- Lines 64-65: Form the transformation matrix; implemented by `A=diag(criss)+fliplr(diag(cross))`.
- Lines 67-68: Transform the coefficients; implemented by `Bkq=transpose(Bkq'*A)`.

### Key state/data transformations

- Lines 45: computes `a{1}` using `a{1}=[1/sqrt(2) 1 1/sqrt(2)]'`.
- Lines 46: computes `a{2}` using `a{2}=[1 1/2 sqrt(6) 1/2 1]'`.
- Lines 47: computes `a{3}` using `a{3}=[sqrt(2) 1/sqrt(3) sqrt(10/3) sqrt(10) sqrt(10/3) 1/sqrt(3) sqrt(2)]'`.
- Lines 48: computes `a{4}` using `a{4}=[2 1/sqrt(2) sqrt(7) sqrt(7/2) 2*sqrt(70) sqrt(7/2) sqrt(7) 1/sqrt(2) 2]'`.
- Lines 49: computes `a{5}` using `a{5}=[2*sqrt(2) 2/sqrt(5) 6*sqrt(2/5) sqrt(3/5) 2*sqrt(21/5) 6*sqrt(14) 2*sqrt(21/5) sqrt(3/5) 6*sqrt(2/5) 2/sqrt(5) 2*sqrt(2)]'`.
- Lines 50: computes `a{6}` using `a{6}=[4 2/sqrt(3) 4*sqrt(11/6) 2*sqrt(11/5) 4*sqrt(11/5) sqrt(22) 4*sqrt(231) sqrt(22) 4*sqrt(11/5) 2*sqrt(11/5) 4*sqrt(11/6) 2/sqrt(3) 4]'`.
- Lines 51: computes `a{7}` using `a{7}=sqrt([32 16/7 416/7 104/7 4576/7 2288/7 13728/7 6864 13728/7 2288/7 4576/7 104/7 416/7 16/7 32]')`.
- Lines 52: computes `a{8}` using `a{8}=sqrt([64 4 120 20/7 1040/7 156/7 1144/7 2860 823680 2860 1144/7 156/7 1040/7 20/7 120 4 64]')`.
- Lines 53: computes `a{9}` using `a{9}=sqrt([128 64/9 2176/9 136/3 2720/9 272/63 7072/21 12771/65 570569/33 1555840 570569/33 12771/65 7072/21 272/63 2720/9 136/3 2176/9 64/9 128]')`.
- Lines 54: computes `a{10}` using `a{10}=sqrt([256 64/5 2432/5 1216/15 82688/15 5168/3 10336/15 5168/15 537472/15 134368/5 11824384 134368/5 537472/15 5168/15 10336/15 5168/3 82688/15 1216/15 2432/5 64/5…`.
- Lines 55: computes `a{11}` using `a{11}=sqrt([512 256/11 10752/11 896/55 68096/55 34048/11 1157632/33 10336/33 82688/55 289408/55 2052166/3 22573824 2052166/3 289408/55 82688/55 10336/33 1157632/33 34048…`.
- Lines 56: computes `a{12}` using `a{12}=sqrt([1024 128/3 5888/3 2944/11 82432/33 20608/33 783104/11 55936/99 1613669/21 475456/11 950912/33 3328192/3 692263936 3328192/3 950912/33 475456/11 1613669/21 55…`.
- Lines 59: computes `criss` using `criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1) ].*a{k}`.
- Lines 62: computes `cross` using `cross=[+1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`.
- Lines 65: computes `A` using `A=diag(criss)+fliplr(diag(cross))`.
- Lines 68: computes `Bkq` using `Bkq=transpose(Bkq'*A)`.

### Local helper functions

- Line 73: `grumble()` — `function grumble(k,Bkq)`. K.W.H. Stevens has done a great disservice to Magnetic Resonance by
  - Representative operation: `if (~isnumeric(k))||(~isreal(k))||(~isfinite(k))|| (~isscalar(k))||(mod(k,1)~=0)||(k<1)||(k>12)`.
  - Representative operation: `(~isscalar(k))||(mod(k,1)~=0)||(k<1)||(k>12)`.

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
- Note: the scaling factors for ranks 7 to 12 were obtained by
- projecting the operators returned by stevens.m onto the
- operators returned by irr_sph_ten.m; the squared factors
- are rational numbers that are the same for every spin
- multiplicity and reproduce the published ranks 1 to 6.

## Implementation structure

- Transforms the coefficients in front of Stevens operators, as
- produced by stevens.m, into the coefficients before the irredu-
- cible spherical tensor operators, as produced by irr_sph_ten.m
- function. Works up to 12th spherical rank. Source for ranks
- up to 6:
- Bkq=stev2sph(k,Bkq)
- k -the spherical rank in question
- Bkq -a column of 2k+1 real coefficients
- in front of Stevens operators, in
- increasing order of projections
- Bkq -a column of 2k+1 complex coefficients
- in front of irreducible spherical

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fliplr()`, `transpose()`, `isscalar()`, `any()`, `iscolumn()`.
