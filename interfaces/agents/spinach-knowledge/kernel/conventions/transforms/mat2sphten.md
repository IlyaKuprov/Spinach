# kernel/conventions/transforms/mat2sphten.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/mat2sphten.m`
- Signature: `[rank0,rank1,rank2]=mat2sphten(M)`
- Total lines: 88

## Purpose

Converts a 3x3 interaction matrix into the irreducible spherical tensor notation: one rank 0 component, three rank 1 components and five rank 2 components to the total of nine independent components. The conventions are matched to Equation (22) of the paper by Len Mueller: The components are listed in the following order: rank 0: (0,0) rank 1: (1,1) (1,0) (1,-1) rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2) and are returne

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Adapt to empty matrices; implemented by `if isempty(M), M=zeros(3); end`.
- Lines 46-47: Check consistency; implemented by `grumble(M)`.
- Lines 49-50: Set result dimensions; implemented by `rank1=zeros(3,1); rank2=zeros(5,1)`.
- Lines 52-53: Rank 0 component; implemented by `rank0=trace(M)/3`.
- Lines 55-56: Rank 1 components; implemented by `rank1(1)=-(1/2)*(M(3,1)-M(1,3)-1i*(M(3,2)-M(2,3)))`.
- Lines 60-61: Rank 2 components; implemented by `rank2(1)=+(1/2)*(M(1,1)-M(2,2)-1i*(M(1,2)+M(2,1)))`.

### Control flow inferred from the code

- Line 44: conditional branch on `isempty(M), M=zeros(3); end`.

### Key state/data transformations

- Lines 50: computes `rank1` using `rank1=zeros(3,1); rank2=zeros(5,1)`.
- Lines 53: computes `rank0` using `rank0=trace(M)/3`.
- Lines 56: computes `rank1(1)` using `rank1(1)=-(1/2)*(M(3,1)-M(1,3)-1i*(M(3,2)-M(2,3)))`.
- Lines 57: computes `rank1(2)` using `rank1(2)=+(1i/sqrt(2))*(M(1,2)-M(2,1))`.
- Lines 58: computes `rank1(3)` using `rank1(3)=-(1/2)*(M(3,1)-M(1,3)+1i*(M(3,2)-M(2,3)))`.
- Lines 61: computes `rank2(1)` using `rank2(1)=+(1/2)*(M(1,1)-M(2,2)-1i*(M(1,2)+M(2,1)))`.
- Lines 62: computes `rank2(2)` using `rank2(2)=-(1/2)*(M(1,3)+M(3,1)-1i*(M(2,3)+M(3,2)))`.
- Lines 63: computes `rank2(3)` using `rank2(3)=+(1/sqrt(6))*(2*M(3,3)-M(1,1)-M(2,2))`.
- Lines 64: computes `rank2(4)` using `rank2(4)=+(1/2)*(M(1,3)+M(3,1)+1i*(M(2,3)+M(3,2)))`.
- Lines 65: computes `rank2(5)` using `rank2(5)=+(1/2)*(M(1,1)-M(2,2)+1i*(M(1,2)+M(2,1)))`.

### Local helper functions

- Line 70: `grumble()` — `function grumble(M)`. When I was 13 I think -Hewlett and Packard were my idols -I called up Bill Hewlett because he lived in Palo Alto and there were no unlisted
  - Representative operation: `if (~isnumeric(M))||(~isreal(M))|| (~ismatrix(M))||any(size(M)~=[3 3])`.
  - Representative operation: `(~ismatrix(M))||any(size(M)~=[3 3])`.

## Syntax

```matlab
[rank0,rank1,rank2]=mat2sphten(M)
```

## Parameters / inputs

- M -3x3 interaction tensor

## Outputs

- rank0 -a single number giving the coefficient of T(0,0) in
- the spherical tensor expansion.
- rank1 -a row vector with three numbers giving the coeffici-
- ents of T(1,1), T(1,0) and T(1,-1) in the spherical
- tensor expansion.
- rank2 -a row vector with five numbers giving the coeffici-
- ents of T(2,2), T(2,1), T(2,0), T(2,-1) and T(2,-2)
- in the spherical tensor expansion.

## Implementation structure

- Converts a 3x3 interaction matrix into the irreducible spherical tensor
- notation: one rank 0 component, three rank 1 components and five rank 2
- components to the total of nine independent components. The conventions
- are matched to Equation (22) of the paper by Len Mueller:
- The components are listed in the following order:
- rank 0: (0,0)
- rank 1: (1,1) (1,0) (1,-1)
- rank 2: (2,2) (2,1) (2,0) (2,-1) (2,-2)
- and are returned as coefficients in front of the corresponding irredu-
- cible spherical tensor operators returned by irr_sph_ten.m function.
- [rank0,rank1,rank2]=mat2sphten(M)
- M -3x3 interaction tensor

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rank1()`, `rank2()`, `ismatrix()`, `any()`.
