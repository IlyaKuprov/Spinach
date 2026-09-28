# kernel/conventions/transforms/anas2mat.m

- Signature: `M=anas2mat(iso,an,as,alp,bet,gam)`

## Purpose

Converts the isotropic part, anisotropy, and asymmetry of a 3x3 interaction tensor into a matrix using Euler angles specified in radians.

## Physical / mathematical content

The function computes principal values `xx`, `yy`, and `zz` from `iso`, `an`, and `as`, then rotates the diagonal tensor into the specified orientation.

## Numerical / algorithmic content

It sets `ra=2*an/3`, `zz=iso+ra`, `yy=iso-ra*(1-as)/2`, and `xx=iso-ra*(1+as)/2`. With `R=euler2dcm(alp,bet,gam)`, it returns `M=R*diag([xx yy zz])*R'`.

## Parameters / inputs

- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- an -interaction anisotropy, defined as zz-(xx+yy)/2
- in terms of eigenvalues
- as -interaction asymmetry, defined as (yy-xx)/(zz-iso)
- in terms of eigenvalues
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians
- gam -gamma Euler angle in radians

## Outputs

- M -3x3 matrix

## Implementation structure

The function checks that all six inputs are real numeric scalars, computes the principal values, and applies the Euler-angle rotation.
