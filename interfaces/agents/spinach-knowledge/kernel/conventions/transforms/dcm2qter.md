# kernel/conventions/transforms/dcm2qter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/dcm2qter.m`
- Signature: `q=dcm2qter(dcm)`
- Total lines: 85

## Purpose

Converts a direction cosine matrix in the active convention of euler2dcm.m function into a unit quaternion. Syntax: q=dcm2qter(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(dcm)`.
- Lines 30-31: Shepperd invariants of the four candidate pivots; implemented by `piv=[1+dcm(1,1)+dcm(2,2)+dcm(3,3)`.
- Lines 36-37: Convert using the best-conditioned pivot; implemented by `[~,best]=max(piv)`.
- Lines 61-62: Resolve the double cover towards a non-negative scalar part; implemented by `if q.u<0`.
- Lines 66-67: Normalise the quaternion; implemented by `qnorm=norm([q.u q.i q.j q.k],2)`.

### Control flow inferred from the code

- Line 38: dispatches on `best`; cases `1`, `2`, `3`, `4`.
- Line 62: conditional branch on `q.u<0`.

### Key state/data transformations

- Lines 31: computes `piv` using `piv=[1+dcm(1,1)+dcm(2,2)+dcm(3,3)`.
- Lines 37: computes `[~,best]` using `[~,best]=max(piv)`.
- Lines 40: computes `q.u` using `q.u=sqrt(piv(1))/2`.
- Lines 41: computes `q.i` using `q.i=(dcm(3,2)-dcm(2,3))/(4*q.u)`.
- Lines 42: computes `q.j` using `q.j=(dcm(1,3)-dcm(3,1))/(4*q.u)`.
- Lines 43: computes `q.k` using `q.k=(dcm(2,1)-dcm(1,2))/(4*q.u)`.
- Lines 67: computes `qnorm` using `qnorm=norm([q.u q.i q.j q.k],2)`.

### Local helper functions

- Line 74: `grumble()` — `function grumble(dcm)`.
  - Representative operation: `if (~isnumeric(dcm))||(~isreal(dcm))||(~isequal(size(dcm),[3 3]))`.
  - Representative operation: `error('dcm must be a real 3x3 matrix.')`.

## Parameters / inputs

- dcm -directional cosine matrix, a 3x3 orthogonal
- matrix with unit determinant

## Outputs

- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion,
- normalised to q.u greater than or equal to zero
- Note: quaternions double-cover rotations; of the two candi-
- dates q and -q this function returns the one with the
- non-negative scalar part.

## Implementation structure

- Converts a direction cosine matrix in the active convention of
- euler2dcm.m function into a unit quaternion. Syntax:
- q=dcm2qter(dcm)
- dcm -directional cosine matrix, a 3x3 orthogonal
- matrix with unit determinant
- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion,
- normalised to q.u greater than or equal to zero
- Note: quaternions double-cover rotations; of the two candi-
- dates q and -q this function returns the one with the
- non-negative scalar part.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dcm()`, `piv()`, `isequal()`.
