# kernel/conventions/transforms/spsk2mat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/spsk2mat.m`
- Signature: `M=spsk2mat(iso,sp,sk,alp,bet,gam)`
- Total lines: 77

## Purpose

Converts span and skew representation of a 3x3 interaction tensor (Herzfeld-Berger convention) into the corresponding matrix. Euler angles should be specified in radians. Syntax: M=spsk2mat(iso,sp,sk,alp,bet,gam)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- sp -interaction span, defined as the difference
- between the largest and the smallest eigenvalue
- sk -interaction skew, defined as 3*(yy-iso)/sp
- where yy is the middle eigenvalue
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians
- gam -gamma Euler angle in radians

## Outputs

- M -3x3 matrix
- Note: the reverse transformation is ill-defined.

## Implementation structure

- Converts span and skew representation of a 3x3 interaction tensor
- (Herzfeld-Berger convention) into the corresponding matrix. Euler
- angles should be specified in radians. Syntax:
- M=spsk2mat(iso,sp,sk,alp,bet,gam)
- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- sp -interaction span, defined as the difference
- between the largest and the smallest eigenvalue
- sk -interaction skew, defined as 3*(yy-iso)/sp
- where yy is the middle eigenvalue
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`, `isscalar()`.
