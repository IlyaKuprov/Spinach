# kernel/indexing/lin2kq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lin2kq.m`
- Signature: `[K,Q]=lin2kq(N,I,idx_base)`
- Total lines: 101

## Purpose

Converts linear serpentine indexing of matrices into their k,q indexing. In base-1 indexing convention: (1,1)(1,2)(1,3) (1)(3)(6) (2,1)(2,2)(2,3) <=> (2)(5)(8) (3,1)(3,2)(3,3) (4)(7)(9) and in base 0 indexing convention: (0,0)(0,1)(0,2) (0)(2)(5) (1,0)(1,1)(1,2) <=> (1)(4)(7) (2,0)(2,1)(2,2) (3)(6)(8)

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(N,I,idx_base)`.
- Lines 44-45: Sepentine matrix; implemented by `S=serpentine(N,idx_base)`.
- Lines 47-48: Direct look-up; implemented by `K=zeros(size(I))`.
- Lines 54-55: 0-base indexing; implemented by `for n=1:numel(I)`.
- Lines 62-63: 1-base indexing; implemented by `for n=1:numel(I)`.
- Lines 69-70: Complain and bomb out; implemented by `error('unsupported indexing base.')`.

### Control flow inferred from the code

- Line 50: dispatches on `idx_base`; cases `0`, `1`.
- Line 55: `for` loop over `n=1:numel(I)`.
- Line 63: `for` loop over `n=1:numel(I)`.

### Key state/data transformations

- Lines 45: computes `S` using `S=serpentine(N,idx_base)`.
- Lines 48: computes `K` using `K=zeros(size(I))`.
- Lines 49: computes `Q` using `Q=zeros(size(I))`.

### Local helper functions

- Line 77: `grumble()` — `function grumble(N,I,idx_base)`.
  - Representative operation: `if (~isnumeric(idx_base))||(~isreal(idx_base))|| (~isscalar(idx_base))||(~ismember(idx_base,[0 1]))`.
  - Representative operation: `(~isscalar(idx_base))||(~ismember(idx_base,[0 1]))`.

## Syntax

```matlab
[K,Q]=lin2kq(N,I,idx_base)
```

## Parameters / inputs

- N -matrix dimension
- I -linear serpentine indices, an
- array of integers
- idx_base -indexing base, 0 or 1

## Outputs

- K -row indices, an array of the
- same size as I
- Q -col indices, an array of the
- same size as I

## Implementation structure

- Converts linear serpentine indexing of matrices into their
- k,q indexing. In base-1 indexing convention:
- (1,1)(1,2)(1,3) (1)(3)(6)
- (2,1)(2,2)(2,3) <=> (2)(5)(8)
- (3,1)(3,2)(3,3) (4)(7)(9)
- and in base 0 indexing convention:
- (0,0)(0,1)(0,2) (0)(2)(5)
- (1,0)(1,1)(1,2) <=> (1)(4)(7)
- (2,0)(2,1)(2,2) (3)(6)(8)
- [K,Q]=lin2kq(N,I,idx_base)
- N -matrix dimension
- I -linear serpentine indices, an

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `serpentine()`, `isscalar()`, `ismember()`, `any()`, `element()`.
