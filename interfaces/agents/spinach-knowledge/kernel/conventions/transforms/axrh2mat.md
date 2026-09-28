# kernel/conventions/transforms/axrh2mat.m

- Signature: `M=axrh2mat(iso,ax,rh,alp,bet,gam)`

## Purpose

Converts axiality and rhombicity representation of a 3x3 interaction tensor into the corresponding matrix. Syntax: M=axrh2mat(iso,ax,rh,alp,bet,gam)

## Physical / mathematical content

The principal values are `xx=iso-(ax+3*rh)/6`, `yy=iso-(ax-3*rh)/6`, and `zz=iso+ax/3`. The Euler-angle rotation maps the diagonal tensor into the requested coordinate frame.
## Numerical / algorithmic content

The function computes the three principal values, obtains a rotation matrix with `euler2dcm(alp,bet,gam)`, and forms `M=R*diag([xx yy zz])*R'`. It then replaces `M` with `(M+M')/2` to enforce symmetry.
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

A local `grumble` function checks that all inputs are real numeric scalars, that `rh` is nonnegative, and that `ax` is at least `rh`. The main function then computes the principal values, rotates the diagonal matrix, and symmetrizes the result.
