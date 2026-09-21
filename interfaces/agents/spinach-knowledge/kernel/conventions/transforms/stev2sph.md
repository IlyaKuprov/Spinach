# kernel/conventions/transforms/stev2sph.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/stev2sph.m`
- Signature: `Bkq=stev2sph(k,Bkq)`
- Total lines: 79

## Purpose

Transforms the coefficients in front of Stevens operators, as produced by stevens.m, into the coefficients before the irredu- cible spherical tensor operators, as produced by irr_sph_ten.m function. Works up to 6th spherical rank. Source:

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Implementation structure

- Transforms the coefficients in front of Stevens operators, as
- produced by stevens.m, into the coefficients before the irredu-
- cible spherical tensor operators, as produced by irr_sph_ten.m
- function. Works up to 6th spherical rank. Source:
- Bkq=stev2sph(k,Bkq)
- k -the spherical rank in question
- Bkq -a column of 2k+1 real coefficients
- in front of Stevens operators, in
- increasing order of projections
- Bkq -a column of 2k+1 complex coefficients
- in front of irreducible spherical
- tensor operators, in decreasing order

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fliplr()`, `transpose()`, `isscalar()`, `any()`, `iscolumn()`.
