# kernel/conventions/transforms/ias2mat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/ias2mat.m`
- Signature: `C=ias2mat(a,d,A)`
- Total lines: 67

## Purpose

Reconstruction of a 3x3 real interaction matrix C between real vectors u and v from its isotropic-antisymmetric-symmetric de- composition: a*(u'*v) + d'*cross(u,v) + u'*A*v = u'*C*v

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(a,d,A)`.
- Lines 32-34: Reconstruct the matrix; implemented by `C=a*eye(3,3)+ [ 0 d(3) -d(2)`.

### Key state/data transformations

- Lines 33-34: computes `C` using `C=a*eye(3,3)+ [ 0 d(3) -d(2)`.

### Local helper functions

- Line 41: `grumble()` — `function grumble(a,d,A)`.
  - Representative operation: `if (~isnumeric(a))||(~isreal(a))||(~isscalar(a))`.
  - Representative operation: `error('a must be a real scalar.')`.

## Syntax

```matlab
C=ias2mat(a,d,A)
```

## Parameters / inputs

- a -scalar component
- d -antisymmetric coupling vector
- A -symmetric coupling matrix

## Outputs

- C -real 3x3 matrix

## Implementation structure

- Reconstruction of a 3x3 real interaction matrix C between real
- vectors u and v from its isotropic-antisymmetric-symmetric de-
- composition:
- a*(u'*v) + d'*cross(u,v) + u'*A*v = u'*C*v
- C=ias2mat(a,d,A)
- a -scalar component
- d -antisymmetric coupling vector
- A -symmetric coupling matrix
- C -real 3x3 matrix
- Check consistency
- Reconstruct the matrix
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
