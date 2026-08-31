# kernel/conventions/transforms/anas2mat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/anas2mat.m`
- Signature: `M=anas2mat(iso,an,as,alp,bet,gam)`
- Total lines: 86

## Purpose

Converts anisotropy and asymmetry representation of a 3x3 interaction tensor (Haeberlen-Mehring convention) into the corresponding matrix. Euler angles should be specified in radians. Syntax: M=anas2mat(iso,an,as,alp,bet,gam)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(iso,an,as,alp,bet,gam)`.
- Lines 37-38: Compute reduced anisotropy; implemented by `ra=2*an/3`.
- Lines 40-41: Compute eigenvalues; implemented by `zz=iso+ra`.
- Lines 45-46: Rotate the matrix; implemented by `R=euler2dcm(alp,bet,gam)`.

### Key state/data transformations

- Lines 38: computes `ra` using `ra=2*an/3`.
- Lines 41: computes `zz` using `zz=iso+ra`.
- Lines 42: computes `yy` using `yy=iso-ra*(1-as)/2`.
- Lines 43: computes `xx` using `xx=iso-ra*(1+as)/2`.
- Lines 46: computes `R` using `R=euler2dcm(alp,bet,gam)`.
- Lines 47: computes `M` using `M=R*diag([xx yy zz])*R'`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(iso,an,as,alp,bet,gam)`. Why is it that beauty is no longer the standard by which we judge things,
  - Representative operation: `if (~isnumeric(iso))||(~isreal(iso))||(~isscalar(iso))|| (~isnumeric(an))||(~isreal(an))||(~isscalar(an))|| (~isnumeric(as))||(~isreal(as))||(~isscalar(as))|| (~isnumeri…`.
  - Representative operation: `(~isnumeric(an))||(~isreal(an))||(~isscalar(an))|| (~isnumeric(as))||(~isreal(as))||(~isscalar(as))|| (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))|| (~isnumeric(b…`.

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

- Converts anisotropy and asymmetry representation of a 3x3 interaction
- tensor (Haeberlen-Mehring convention) into the corresponding matrix.
- Euler angles should be specified in radians. Syntax:
- M=anas2mat(iso,an,as,alp,bet,gam)
- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- an -interaction anisotropy, defined as zz-(xx+yy)/2
- in terms of eigenvalues
- as -interaction asymmetry, defined as (yy-xx)/(zz-iso)
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians
- gam -gamma Euler angle in radians

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`, `isscalar()`.
