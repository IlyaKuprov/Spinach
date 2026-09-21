# kernel/conventions/transforms/stev2sph.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/stev2sph.m`
- Signature: `Bkq=stev2sph(k,Bkq)`
- Total lines: 98

## Purpose

Transforms the coefficients in front of Stevens operators, as produced by stevens.m, into the coefficients before the irredu- cible spherical tensor operators, as produced by irr_sph_ten.m function. Works up to 12th spherical rank. Source for ranks up to 6:

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(k,Bkq)`.
- Lines 49-50: Catalog the stupid scaling factors; implemented by `a{1}=[1/sqrt(2) 1 1/sqrt(2)]'`.
- Lines 64-65: Form the transformation matrix diagonal; implemented by `criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1) ].*a{k}`.
- Lines 67-68: Form the transformation matrix antidiagonal; implemented by `cross=[+1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`.
- Lines 70-71: Form the transformation matrix; implemented by `A=diag(criss)+fliplr(diag(cross))`.
- Lines 73-74: Transform the coefficients; implemented by `Bkq=transpose(Bkq'*A)`.

### Key state/data transformations

- Lines 50: computes `a{1}` using `a{1}=[1/sqrt(2) 1 1/sqrt(2)]'`.
- Lines 51: computes `a{2}` using `a{2}=[1 1/2 sqrt(6) 1/2 1]'`.
- Lines 52: computes `a{3}` using `a{3}=[sqrt(2) 1/sqrt(3) sqrt(10/3) sqrt(10) sqrt(10/3) 1/sqrt(3) sqrt(2)]'`.
- Lines 53: computes `a{4}` using `a{4}=[2 1/sqrt(2) sqrt(7) sqrt(7/2) 2*sqrt(70) sqrt(7/2) sqrt(7) 1/sqrt(2) 2]'`.
- Lines 54: computes `a{5}` using `a{5}=[2*sqrt(2) 2/sqrt(5) 6*sqrt(2/5) sqrt(3/5) 2*sqrt(21/5) 6*sqrt(14) 2*sqrt(21/5) sqrt(3/5) 6*sqrt(2/5) 2/sqrt(5) 2*sqrt(2)]'`.
- Lines 55: computes `a{6}` using `a{6}=[4 2/sqrt(3) 4*sqrt(11/6) 2*sqrt(11/5) 4*sqrt(11/5) sqrt(22) 4*sqrt(231) sqrt(22) 4*sqrt(11/5) 2*sqrt(11/5) 4*sqrt(11/6) 2/sqrt(3) 4]'`.
- Lines 56: computes `a{7}` using `a{7}=sqrt([32 16/7 416/7 104/7 4576/7 2288/7 13728/7 6864 13728/7 2288/7 4576/7 104/7 416/7 16/7 32]')`.
- Lines 57: computes `a{8}` using `a{8}=sqrt([64 4 120 20/7 1040/7 156/7 1144/7 2860 823680 2860 1144/7 156/7 1040/7 20/7 120 4 64]')`.
- Lines 58-59: computes `a{9}` using `a{9}=sqrt([128 64/9 2176/9 136/3 2720/9 272/63 7072/21 28742418432/146289025 2529332822016/146289025 1555840 2529332822016/146289025 28742418432/146289025 7072/21 272/63…`.
- Lines 60: computes `a{10}` using `a{10}=sqrt([256 64/5 2432/5 1216/15 82688/15 5168/3 10336/15 5168/15 537472/15 134368/5 11824384 134368/5 537472/15 5168/15 10336/15 5168/3 82688/15 1216/15 2432/5 64/5…`.
- Lines 61: computes `a{11}` using `a{11}=sqrt([512 256/11 10752/11 896/55 68096/55 34048/11 1157632/33 10336/33 82688/55 289408/55 7524608/11 22573824 7524608/11 289408/55 82688/55 10336/33 1157632/33 340…`.
- Lines 62: computes `a{12}` using `a{12}=sqrt([1024 128/3 5888/3 2944/11 82432/33 20608/33 783104/11 55936/99 7607296/99 475456/11 950912/33 3328192/3 692263936 3328192/3 950912/33 475456/11 7607296/99 55…`.
- Lines 65: computes `criss` using `criss=[-1i*(-1).^(k:-1:1)'; 1; ones(k,1) ].*a{k}`.
- Lines 68: computes `cross` using `cross=[+1i*ones(k,1); 0; (-1).^(1:k)'].*a{k}`.
- Lines 71: computes `A` using `A=diag(criss)+fliplr(diag(cross))`.
- Lines 74: computes `Bkq` using `Bkq=transpose(Bkq'*A)`.

### Local helper functions

- Line 79: `grumble()` — `function grumble(k,Bkq)`. K.W.H. Stevens has done a great disservice to Magnetic Resonance by
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
- Note: the squared scaling factors for ranks 7 to 12 are exact
- rationals computed from the integer coefficient table of
- stevens.m (Ryabov, J. Magn. Reson. 140, 141 (1999)) and
- the normalisation of irr_sph_ten.m: 2^(k-2)*P(k,q)/C(k,q)^2
- for q>0 and 2^k*P(k,0)/C(k,0)^2 for q=0, where P(k,q) is
- the product of (k+p)(k-p+1) over p from q+1 to k, C(k,q)
- is the stevens.m coefficient (halved for even k and odd q),
- and the same expression reproduces the published ranks 1
- to 6. The Ryabov table has a cluster of large primes at
- rank 9 projections 1 and 2, hence the denominators there.

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
