# kernel/conventions/transforms/axrh2mat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/axrh2mat.m`
- Signature: `M=axrh2mat(iso,ax,rh,alp,bet,gam)`
- Total lines: 72

## Purpose

Converts axiality and rhombicity representation of a 3x3 interaction tensor into the corresponding matrix. Syntax: M=axrh2mat(iso,ax,rh,alp,bet,gam)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(iso,ax,rh,alp,bet,gam)`.
- Lines 40-41: Compute eigenvalues; implemented by `xx=iso-(ax+3*rh)/6`.
- Lines 45-46: Rotate the interaction; implemented by `R=euler2dcm(alp,bet,gam)`.
- Lines 49-50: Tidy up; implemented by `M=(M+M')/2`.

### Key state/data transformations

- Lines 41: computes `xx` using `xx=iso-(ax+3*rh)/6`.
- Lines 42: computes `yy` using `yy=iso-(ax-3*rh)/6`.
- Lines 43: computes `zz` using `zz=iso+ax/3`.
- Lines 46: computes `R` using `R=euler2dcm(alp,bet,gam)`.
- Lines 47: computes `M` using `M=R*diag([xx yy zz])*R'`.

### Local helper functions

- Line 55: `grumble()` — `function grumble(iso,ax,rh,alp,bet,gam)`.
  - Representative operation: `if (~isnumeric(iso))||(~isreal(iso))||(~isscalar(iso))|| (~isnumeric(ax))||(~isreal(ax))||(~isscalar(ax))|| (~isnumeric(rh))||(~isreal(rh))||(~isscalar(rh))|| (~isnumeri…`.
  - Representative operation: `(~isnumeric(ax))||(~isreal(ax))||(~isscalar(ax))|| (~isnumeric(rh))||(~isreal(rh))||(~isscalar(rh))|| (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))|| (~isnumeric(b…`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`, `isscalar()`.
