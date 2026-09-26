# kernel/conventions/transforms/axrh2mat.m

- Signature: `M=axrh2mat(iso,ax,rh,alp,bet,gam)`

## Purpose

Converts axiality and rhombicity representation of a 3x3 interaction tensor into the corresponding matrix. Syntax: M=axrh2mat(iso,ax,rh,alp,bet,gam)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- ax -interaction axiality, defined as 2*zz-(xx+yy)
- in terms of eigenvalues (Mehring order, that
- is xx<=yy<=zz)
- rh -interaction rhombicity, defined as (yy-xx) in
- terms of eigenvalues (Mehring order, that is
- xx<=yy<=zz)
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians
- gam -gamma Euler angle in radians

## Outputs

- M -3x3 matrix
- Note: the inverse transformation is ill-defined.

## Implementation structure

- Converts axiality and rhombicity representation of a 3x3 interaction
- tensor into the corresponding matrix. Syntax:
- M=axrh2mat(iso,ax,rh,alp,bet,gam)
- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- ax -interaction axiality, defined as 2*zz-(xx+yy)
- in terms of eigenvalues (Mehring order, that
- is xx<=yy<=zz)
- rh -interaction rhombicity, defined as (yy-xx) in
- terms of eigenvalues (Mehring order, that is
- xx<=yy<=zz)
- alp -alpha Euler angle in radians
