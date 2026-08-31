# kernel/indexing/kq2lin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/kq2lin.m`
- Signature: `I=kq2lin(N,K,Q,idx_base)`
- Total lines: 104

## Purpose

Converts k,q indexing of matrices into their linear serpentine indexing. In base 1 indexing convention: (1,1)(1,2)(1,3) (1)(3)(6) (2,1)(2,2)(2,3) <=> (2)(5)(8) (3,1)(3,2)(3,3) (4)(7)(9) and in base 0 indexing convention: (0,0)(0,1)(0,2) (0)(2)(5) (1,0)(1,1)(1,2) <=> (1)(4)(7) (2,0)(2,1)(2,2) (3)(6)(8)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(N,K,Q,idx_base)`.
- Lines 44-45: Sepentine matrix; implemented by `S=serpentine(N,idx_base)`.
- Lines 47-48: Direct look-up; implemented by `I=zeros(size(K))`.
- Lines 53-54: 0-base indexing; implemented by `for n=1:numel(K)`.
- Lines 60-61: 1-base indexing; implemented by `for n=1:numel(K)`.
- Lines 67-68: Complain and bomb out; implemented by `error('unsupported indexing base.')`.

### Control flow inferred from the code

- Line 49: dispatches on `idx_base`; cases `0`, `1`.
- Line 54: `for` loop over `n=1:numel(K)`.
- Line 61: `for` loop over `n=1:numel(K)`.

### Key state/data transformations

- Lines 45: computes `S` using `S=serpentine(N,idx_base)`.
- Lines 48: computes `I` using `I=zeros(size(K))`.
- Lines 55: computes `I(n)` using `I(n)=S(K(n)+1,Q(n)+1)`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(N,K,Q,idx_base)`.
  - Representative operation: `if (~isnumeric(idx_base))||(~isreal(idx_base))|| (~isscalar(idx_base))||(~ismember(idx_base,[0 1]))`.
  - Representative operation: `(~isscalar(idx_base))||(~ismember(idx_base,[0 1]))`.

## Syntax

```matlab
I=kq2lin(N,K,Q,idx_base)
```

## Parameters / inputs

- N -matrix dimension, a scalar
- K -first index, positive integer
- array of any size
- Q -second index, positive integer
- array of the same size as K
- idx_base -indexing base, 0 or 1

## Outputs

- I -linear serpentine index, an ar-
- ray of the same size as inputs

## Implementation structure

- Converts k,q indexing of matrices into their linear
- serpentine indexing. In base 1 indexing convention:
- (1,1)(1,2)(1,3) (1)(3)(6)
- (2,1)(2,2)(2,3) <=> (2)(5)(8)
- (3,1)(3,2)(3,3) (4)(7)(9)
- and in base 0 indexing convention:
- (0,0)(0,1)(0,2) (0)(2)(5)
- (1,0)(1,1)(1,2) <=> (1)(4)(7)
- (2,0)(2,1)(2,2) (3)(6)(8)
- I=kq2lin(N,K,Q,idx_base)
- N -matrix dimension, a scalar
- K -first index, positive integer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `serpentine()`, `isscalar()`, `ismember()`, `any()`, `element()`.
