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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(D,E,alp,bet,gam)`.
- Lines 30-31: Compute the matrix in the eigenframe; implemented by `M=[-D/3+E, 0, 0; 0, -D/3-E, 0; 0, 0, 2*D/3]`.
- Lines 33-34: Rotate the molecule; implemented by `R=euler2dcm(alp,bet,gam); M=R*M*R'`.
- Lines 36-37: Tidy up the double precision; implemented by `M=M-eye(3)*trace(M)/3; M=(M+M')/2`.

### Key state/data transformations

- Lines 31: computes `M` using `M=[-D/3+E, 0, 0; 0, -D/3-E, 0; 0, 0, 2*D/3]`.
- Lines 34: computes `R` using `R=euler2dcm(alp,bet,gam); M=R*M*R'`.

### Local helper functions

- Line 42: `grumble()` — `function grumble(D,E,alp,bet,gam)`. To watch the courageous Afghan freedom fighters battle modern arsenals with simple hand-held weapons is an inspiration to
  - Representative operation: `if (~isnumeric(D))||(~isreal(D))||(~isscalar(D))|| (~isnumeric(E))||(~isreal(E))||(~isscalar(E))|| (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))|| (~isnumeric(bet)…`.
  - Representative operation: `(~isnumeric(E))||(~isreal(E))||(~isscalar(E))|| (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))|| (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))|| (~isnumeric(g…`.

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
