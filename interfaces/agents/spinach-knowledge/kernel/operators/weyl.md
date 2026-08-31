# kernel/operators/weyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/weyl.m`
- Signature: `A=weyl(nlevels)`
- Total lines: 78

## Purpose

Weyl boson operators (sparse, see below for normalisa- tion convention) for a bosonic mode with a user-speci- fied population number truncation. Syntax: A=weyl(nlevels)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Ensure internal consistency; implemented by `grumble(nlevels); nlevels=double(nlevels)`.
- Lines 43-44: Creation operator; implemented by `diags=sqrt(1:nlevels)`.
- Lines 48-49: Number operator; implemented by `diags=0:(nlevels-1)`.
- Lines 53-54: Annihilation operator; implemented by `diags=sqrt(0:(nlevels-1))`.
- Lines 58-59: Unit operator; implemented by `A.u=complex(speye(nlevels))`.

### Key state/data transformations

- Lines 41: computes `grumble(nlevels); nlevels` using `grumble(nlevels); nlevels=double(nlevels)`.
- Lines 44: computes `diags` using `diags=sqrt(1:nlevels)`.
- Lines 45: computes `A.c` using `A.c=spdiags(diags',-1,nlevels,nlevels)`.
- Lines 50: computes `A.n` using `A.n=spdiags(diags',0,nlevels,nlevels)`.
- Lines 55: computes `A.a` using `A.a=spdiags(diags',+1,nlevels,nlevels)`.
- Lines 59: computes `A.u` using `A.u=complex(speye(nlevels))`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(nlevels)`. I happen to be a physicist who started life as a mathematician. As a working physicist, I am acutely
  - Representative operation: `if (~isnumeric(nlevels))||(~isreal(nlevels))|| (~isscalar(nlevels))||(mod(nlevels,1)~=0)||(nlevels<1)`.
  - Representative operation: `(~isscalar(nlevels))||(mod(nlevels,1)~=0)||(nlevels<1)`.

## Parameters / inputs

- nlevels -an integer specifying the
- number of population levels

## Outputs

- A.u -unit operator
- A.c -creation operator
- A.a -annihilation operator
- A.n -population number operator
- Note: the matrices are normalised to obey the following
- relations for all energy level counts
- A.c*A.a=A.n, [A.n,A.c]=A.c
- [A.n,A.a]=-A.a, [A.a,A.c]=A.u
- except for the edge state at which [A.a,A.c] ele-
- ment is (1-nlevels), this is unavoidable.
- Note: arrays are declared complex at build time to avoid
- expensive reallocation operations later on.

## Implementation structure

- Weyl boson operators (sparse, see below for normalisa-
- tion convention) for a bosonic mode with a user-speci-
- fied population number truncation. Syntax:
- A=weyl(nlevels)
- nlevels -an integer specifying the
- number of population levels
- A.u -unit operator
- A.c -creation operator
- A.a -annihilation operator
- A.n -population number operator
- Note: the matrices are normalised to obey the following
- relations for all energy level counts

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `double()`, `spdiags()`, `complex()`, `speye()`, `isscalar()`.
