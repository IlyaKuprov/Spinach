# kernel/conventions/transforms/zfs2mat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/zfs2mat.m`
- Signature: `M=zfs2mat(D,E,alp,bet,gam)`
- Total lines: 61

## Purpose

Converts D and E zero-field splitting parameters described in the abstract of (http://dx.doi.org/10.1063/1.1682294) into a spin interaction matrix. Syntax: M=zfs2mat(D,E,alp,bet,gam)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- D,E -real scalar parameters, Hz
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians
- gam -gamma Euler angle in radians

## Outputs

- M -symmetric 3x3 matrix, Hz

## Implementation structure

- Converts D and E zero-field splitting parameters described in
- the abstract of (http://dx.doi.org/10.1063/1.1682294) into a
- spin interaction matrix. Syntax:
- M=zfs2mat(D,E,alp,bet,gam)
- D,E -real scalar parameters, Hz
- alp -alpha Euler angle in radians
- bet -beta Euler angle in radians
- gam -gamma Euler angle in radians
- M -symmetric 3x3 matrix, Hz
- Check consistency
- Compute the matrix in the eigenframe
- Rotate the molecule

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`, `isscalar()`.
