# kernel/conventions/transforms/mat2ias.m

- Signature: `[a,d,A]=mat2ias(C)`

## Purpose

Isotropic-antisymmetric-symmetric decomposition of a 3x3 real interaction matrix between real vectors u and v: u'*C*v = a*(u'*v) + d'*cross(u,v) + u'*A*v

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

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
