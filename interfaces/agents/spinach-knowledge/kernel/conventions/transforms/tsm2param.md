# kernel/conventions/transforms/tsm2param.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/tsm2param.m`
- Signature: `[ax,rh,angles]=tsm2param(M)`
- Total lines: 79

## Purpose

Attempts to convert a traceless symmetric 3x3 interaction matrix into axiality, rhombicity and three Euler angles. The transformation is un- stable and should be avoided if at all possible: it is always best to just publish the 3x3 matrix as recommended by IUPAC. Syntax: [ax,rh,angles]=tsm2param(M)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- M -3x3 matrix or its five independent elements in the
- order of [Mxx, Mxy, Mxz, Myy, Myz]

## Outputs

- ax -axiality, Mehring order of eigenvalues
- rh -rhombicity, Mehring order of eigenvalues
- angles -Euler angles (one of the eight equivalent
- sets), radians
- Note: Mehring convention has Z as the largest eigenvalue, and X as
- the smallest eigenvalue, this includes signs.

## Implementation structure

- Attempts to convert a traceless symmetric 3x3 interaction matrix into
- axiality, rhombicity and three Euler angles. The transformation is un-
- stable and should be avoided if at all possible: it is always best to
- just publish the 3x3 matrix as recommended by IUPAC. Syntax:
- [ax,rh,angles]=tsm2param(M)
- M - 3x3 matrix or its five independent elements in the
- order of [Mxx, Mxy, Mxz, Myy, Myz]
- ax - axiality, Mehring order of eigenvalues
- rh - rhombicity, Mehring order of eigenvalues
- angles - Euler angles (one of the eight equivalent
- sets), radians
- Note: Mehring convention has Z as the largest eigenvalue, and X as

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `setdiff()`, `dcm2euler()`, `issymmetric()`.
