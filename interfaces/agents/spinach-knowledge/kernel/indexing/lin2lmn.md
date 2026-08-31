# kernel/indexing/lin2lmn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/indexing/lin2lmn.m`
- Signature: `[L,M,N]=lin2lmn(I)`
- Total lines: 66

## Purpose

Converts linear indices of Wigner D functions into L,M,N indices. In the linear indexing convention, Wigner D functions are listed in the order of increasing L rank. Within each L, the functions are listed in the order of decreasing left index M, and, for each M, in the or- der of decreasing N index. One base counting is used: I=1 -> (L=0,M=0,N=0) I=2 -> (L=1,M=1,N=1) I=3 -> (L=1,M=1,N=0), et cetera...

## Physical / mathematical content

- Indexing utilities. These files build and transform compact index maps for basis states, matrix elements, trajectories, and tensor-product structures.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(I)`.
- Lines 37-38: Solve the cubic equation for L; implemented by `big_root=(27*I+sqrt(729*I.^2-3)).^(1/3)`.
- Lines 41-42: Get the left index; implemented by `rank_page_position=I-(4*L.^3-L)/3-1`.
- Lines 45-46: Get the right index; implemented by `N=L+(2*L+1).*(L-M)-rank_page_position`.
- Lines 48-49: Make sure the conversion is correct; implemented by `if nnz(lmn2lin(L,M,N)~=I)>0`.

### Control flow inferred from the code

- Line 49: conditional branch on `nnz(lmn2lin(L,M,N)~=I)>0`.

### Key state/data transformations

- Lines 38: computes `big_root` using `big_root=(27*I+sqrt(729*I.^2-3)).^(1/3)`.
- Lines 39: computes `L` using `L=ceil((3^(1/3)+big_root.^2)./(2*(3^(2/3))*big_root)-1)`.
- Lines 42: computes `rank_page_position` using `rank_page_position=I-(4*L.^3-L)/3-1`.
- Lines 43: computes `M` using `M=L-fix(rank_page_position./(2*L+1))`.
- Lines 46: computes `N` using `N=L+(2*L+1).*(L-M)-rank_page_position`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(I)`. Who gave you the authority to decide which colour Spinach logo was going to be?
  - Representative operation: `if (~isnumeric(I))||(~isreal(I))||any(mod(I(:),1)~=0)||any(I(:)<1)`.
  - Representative operation: `error('all elements of the input array must be positive integers.')`.

## Syntax

```matlab
[L,M,N]=lin2lmn(I)
```

## Parameters / inputs

- I -linear indices of Wigner D functions, with
- I=1 corresponding to L=0, M=0, N=0.

## Outputs

- L -ranks of Wigner D functions
- M -row indices of Wigner D functions
- N -column indices of Wigner D functions

## Implementation structure

- Converts linear indices of Wigner D functions into L,M,N indices. In
- the linear indexing convention, Wigner D functions are listed in the
- order of increasing L rank. Within each L, the functions are listed
- in the order of decreasing left index M, and, for each M, in the or-
- der of decreasing N index. One base counting is used:
- I=1 -> (L=0,M=0,N=0)
- I=2 -> (L=1,M=1,N=1)
- I=3 -> (L=1,M=1,N=0), et cetera...
- [L,M,N]=lin2lmn(I)
- I -linear indices of Wigner D functions, with
- I=1 corresponding to L=0, M=0, N=0.
- L -ranks of Wigner D functions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fix()`, `nnz()`, `lmn2lin()`, `any()`.
