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
