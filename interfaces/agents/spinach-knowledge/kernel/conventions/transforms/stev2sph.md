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
