# kernel/operators/boson_mono.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/boson_mono.m`
- Signature: `B=boson_mono(nlevels)`
- Total lines: 70

## Purpose

Bosonic monomial operators of the following structure: B(k,q)=(Cr^k)*(An^q) obeying the following commutation relations with the po- pulation number operator N: [N,B(k,q)]=(k-q)*B(k,q)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(nlevels)`.
- Lines 37-38: Get generators; implemented by `A=weyl(nlevels)`.
- Lines 40-41: Generate monomials; implemented by `B=cell(nlevels,nlevels)`.
- Lines 48-49: Build the serpentine index; implemented by `[rows,cols]=ndgrid(1:nlevels)`.
- Lines 52-53: Arrange; implemented by `B=B(idx)`.

### Control flow inferred from the code

- Line 42: `for` loop over `k=0:(nlevels-1)`.
- Line 43: `for` loop over `q=0:(nlevels-1)`.

### Key state/data transformations

- Lines 38: computes `A` using `A=weyl(nlevels)`.
- Lines 41: computes `B` using `B=cell(nlevels,nlevels)`.
- Lines 44: computes `B{k+1,q+1}` using `B{k+1,q+1}=(A.c^k)*(A.a^q)`.
- Lines 49: computes `[rows,cols]` using `[rows,cols]=ndgrid(1:nlevels)`.
- Lines 50: computes `[~,idx]` using `[~,idx]=sortrows([rows(:)+cols(:), -rows(:)])`.

### Local helper functions

- Line 58: `grumble()` — `function grumble(nlevels)`. "Oh, so you really exist. I thought Littlewood was a pseudonym Hardy used for his less impor-
  - Representative operation: `if (~isnumeric(nlevels))||(~isreal(nlevels))|| (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.
  - Representative operation: `(~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)`.

## Syntax

```matlab
B=boson_mono(nlevels)
```

## Parameters / inputs

- nlevels -number of bosonic ladder population
- levels, k and q go from 0 to nlevels-1

## Outputs

- B -a cell array with the following numbering
- map between (k,q) and a single index:
- (0,0)(0,1)(0,2) (1)(3)(6)
- (1,0)(1,1)(1,2) <=> (2)(5)(8)
- (2,0)(2,1)(2,2) (4)(7)(9)

## Implementation structure

- Bosonic monomial operators of the following structure:
- B(k,q)=(Cr^k)*(An^q)
- obeying the following commutation relations with the po-
- pulation number operator N:
- [N,B(k,q)]=(k-q)*B(k,q)
- B=boson_mono(nlevels)
- nlevels -number of bosonic ladder population
- levels, k and q go from 0 to nlevels-1
- B -a cell array with the following numbering
- map between (k,q) and a single index:
- (0,0)(0,1)(0,2) (1)(3)(6)
- (1,0)(1,1)(1,2) <=> (2)(5)(8)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `weyl()`, `sortrows()`, `rows()`, `cols()`, `isscalar()`.
