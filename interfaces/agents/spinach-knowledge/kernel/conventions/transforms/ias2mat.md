# kernel/conventions/transforms/ias2mat.m

- Signature: `C=ias2mat(a,d,A)`

## Purpose

Reconstruction of a 3x3 real interaction matrix C between real vectors u and v from its isotropic-antisymmetric-symmetric de- composition: a*(u'*v) + d'*cross(u,v) + u'*A*v = u'*C*v

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

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
