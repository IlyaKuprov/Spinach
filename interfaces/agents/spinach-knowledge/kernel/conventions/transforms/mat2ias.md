# kernel/conventions/transforms/mat2ias.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/mat2ias.m`
- Signature: `[a,d,A]=mat2ias(C)`
- Total lines: 55

## Purpose

Isotropic-antisymmetric-symmetric decomposition of a 3x3 real interaction matrix between real vectors u and v: u'*C*v = a*(u'*v) + d'*cross(u,v) + u'*A*v

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(C)`.
- Lines 31-32: Isotropic part; implemented by `a=trace(C)/3`.
- Lines 34-35: Antisymmetric part; implemented by `d=[C(2,3)-C(3,2)`.
- Lines 39-40: Traceless symmetric part; implemented by `A=(C+C')/2-a*eye(3,3)`.

### Key state/data transformations

- Lines 32: computes `a` using `a=trace(C)/3`.
- Lines 35: computes `d` using `d=[C(2,3)-C(3,2)`.
- Lines 40: computes `A` using `A=(C+C')/2-a*eye(3,3)`.

### Local helper functions

- Line 45: `grumble()` — `function grumble(C)`. Mathematics is the only true metaphysics. Lord Kelvin
  - Representative operation: `if (~isnumeric(C))||(~isreal(C))|| (size(C,1)~=3)||(size(C,2)~=3)`.
  - Representative operation: `(size(C,1)~=3)||(size(C,2)~=3)`.

## Syntax

```matlab
[a,d,A]=mat2ias(C)
```

## Parameters / inputs

- C -real 3x3 matrix

## Outputs

- a -scalar component
- d -antisymmetric coupling vector
- A -symmetric coupling matrix

## Implementation structure

- Isotropic-antisymmetric-symmetric decomposition of a 3x3
- real interaction matrix between real vectors u and v:
- u'*C*v = a*(u'*v) + d'*cross(u,v) + u'*A*v
- [a,d,A]=mat2ias(C)
- C -real 3x3 matrix
- a -scalar component
- d -antisymmetric coupling vector
- A -symmetric coupling matrix
- Check consistency
- Isotropic part
- Antisymmetric part
- Traceless symmetric part

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
