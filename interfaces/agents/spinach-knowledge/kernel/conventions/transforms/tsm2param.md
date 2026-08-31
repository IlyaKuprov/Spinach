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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(M)`.
- Lines 35-36: Assemble the matrix if needed; implemented by `if numel(M)==5`.
- Lines 42-43: Diagonalise the matrix; implemented by `[V,D]=eig(M); D=diag(D)`.
- Lines 45-46: Mehring order eigenvalues; implemented by `[~,IZ]=max(D); [~,IX]=min(D)`.
- Lines 49-50: Compute the invariants; implemented by `ax=2*D(IZ)-(D(IX)+D(IY))`.
- Lines 53-54: Compute Euler angles; implemented by `V=V(:,[IX IY IZ])`.

### Control flow inferred from the code

- Line 36: conditional branch on `numel(M)==5`.

### Key state/data transformations

- Lines 37: computes `M` using `M=[M(1) M(2) M(3)`.
- Lines 43: computes `[V,D]` using `[V,D]=eig(M); D=diag(D)`.
- Lines 46: computes `[~,IZ]` using `[~,IZ]=max(D); [~,IX]=min(D)`.
- Lines 47: computes `IY` using `IY=setdiff([1 2 3],[IZ IX])`.
- Lines 50: computes `ax` using `ax=2*D(IZ)-(D(IX)+D(IY))`.
- Lines 51: computes `rh` using `rh=D(IY)-D(IX)`.
- Lines 54: computes `V` using `V=V(:,[IX IY IZ])`.
- Lines 55: computes `angles` using `angles=dcm2euler(V*det(V))`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(M)`.
  - Representative operation: `if (~isnumeric(M))||(~isreal(M))`.
  - Representative operation: `error('M must be a real numeric array.')`.

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
