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
