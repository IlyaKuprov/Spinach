# kernel/indexing/lmn2lin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lmn2lin.m`
- Signature: `I=lmn2lin(L,M,N)`
- Total lines: 76

## Purpose

Converts L,M,N indices of Wigner D functions into linear indices. In the linear indexing convention, Wigner D functions are listed in the order of increasing L rank. Within each L, the functions are listed in the order of decreasing left index M, and, for each M, in the or- der of decreasing N index. One base counting is used: (L=0,M=0,N=0) -> I=1 (L=1,M=1,N=1) -> I=2 (L=1,M=1,N=0) -> I=3, et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(L,M,N)`.
- Lines 37-38: Get the linear index; implemented by `I=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1`.

### Key state/data transformations

- Lines 38: computes `I` using `I=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1`.

### Local helper functions

- Line 43: `grumble()` — `function grumble(L,M,N)`.
  - Representative operation: `if (~isnumeric(L))||(~isreal(L))||any(mod(L(:),1)~=0)|| (~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)|| (~isnumeric(N))||(~isreal(N))||any(mod(N(:),1)~=0)`.
  - Representative operation: `(~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)|| (~isnumeric(N))||(~isreal(N))||any(mod(N(:),1)~=0)`.

## Syntax

```matlab
I=lmn2lin(L,M,N)
```

## Parameters / inputs

- L -ranks of Wigner D functions
- M -row indices of Wigner D functions
- N -column indices of Wigner D functions

## Outputs

- I -linear indices of Wigner D functions, with
- I=1 corresponding to L=0, M=0, N=0.

## Implementation structure

- Converts L,M,N indices of Wigner D functions into linear indices. In
- the linear indexing convention, Wigner D functions are listed in the
- order of increasing L rank. Within each L, the functions are listed
- in the order of decreasing left index M, and, for each M, in the or-
- der of decreasing N index. One base counting is used:
- (L=0,M=0,N=0) -> I=1
- (L=1,M=1,N=1) -> I=2
- (L=1,M=1,N=0) -> I=3, et cetera...
- I=lmn2lin(L,M,N)
- L -ranks of Wigner D functions
- M -row indices of Wigner D functions
- N -column indices of Wigner D functions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`.
