# kernel/conventions/transforms/qter2euler.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/qter2euler.m`
- Signature: `[alpha,beta,gamma]=qter2euler(q)`
- Total lines: 65

## Purpose

Converts a unit quaternion in the active convention into Euler angles (ZYZ active convention), matching euler2dcm.m function.

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(q)`.
- Lines 33-34: Normalise the quaternions; implemented by `qnorm=sqrt(q.u.^2+q.i.^2+q.j.^2+q.k.^2)`.
- Lines 38-39: Compute the angle sums and differences; implemented by `sum_ag=2*atan2(q.k,q.u); dif_ga=2*atan2(q.i,q.j)`.
- Lines 41-42: Compute the Euler angles; implemented by `beta=2*atan2(sqrt(q.i.^2+q.j.^2),sqrt(q.u.^2+q.k.^2))`.

### Key state/data transformations

- Lines 34: computes `qnorm` using `qnorm=sqrt(q.u.^2+q.i.^2+q.j.^2+q.k.^2)`.
- Lines 35: computes `q.u` using `q.u=q.u./qnorm; q.i=q.i./qnorm`.
- Lines 36: computes `q.j` using `q.j=q.j./qnorm; q.k=q.k./qnorm`.
- Lines 39: computes `sum_ag` using `sum_ag=2*atan2(q.k,q.u); dif_ga=2*atan2(q.i,q.j)`.
- Lines 42: computes `beta` using `beta=2*atan2(sqrt(q.i.^2+q.j.^2),sqrt(q.u.^2+q.k.^2))`.
- Lines 43: computes `alpha` using `alpha=(sum_ag-dif_ga)/2; gamma=(sum_ag+dif_ga)/2`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(q)`.
  - Representative operation: `if ~all(isfield(q,{'i','j','k','u'}))`.
  - Representative operation: `error('quaternion data structure must contain u, i, j, and k fields.')`.

## Syntax

```matlab
[alpha,beta,gamma]=qter2euler(q)
```

## Parameters / inputs

- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; each field may
- be a column vector, in which case the conversion is
- performed elementwise

## Outputs

- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), same shape as the qua-
- ternion component fields
- Note: Euler angles are not unique; the angles returned sa-
- tisfy euler2dcm(alpha,beta,gamma)=qter2dcm(q) with
- beta in the [0,pi] interval.

## Implementation structure

- Converts a unit quaternion in the active convention into Euler
- angles (ZYZ active convention), matching euler2dcm.m function.
- [alpha,beta,gamma]=qter2euler(q)
- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; each field may
- be a column vector, in which case the conversion is
- performed elementwise
- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), same shape as the qua-
- ternion component fields
- Note: Euler angles are not unique; the angles returned sa-
- tisfy euler2dcm(alpha,beta,gamma)=qter2dcm(q) with

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `atan2()`, `all()`, `isfield()`, `iscolumn()`, `any()`, `eps()`.
