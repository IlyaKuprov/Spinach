# kernel/conventions/transforms/mat2axrh.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/mat2axrh.m`
- Signature: `[iso,ax,rh,eigvals]=mat2axrh(M)`
- Total lines: 63

## Purpose

Computes axiality and rhombicity of a symmetric 3x3 interaction tensor from the corresponding matrix. Syntax: [iso,ax,rh]=mat2axrh(M)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(M)`.
- Lines 38-39: Get the eigenvalues; implemented by `[~,eigvals]=eig(M,'vector')`.
- Lines 41-42: Sort in Mehring order; implemented by `eigvals=sort(eigvals,'ascend')`.
- Lines 44-45: Get the interaction parameters; implemented by `ax=2*eigvals(3)-(eigvals(1)+eigvals(2))`.

### Key state/data transformations

- Lines 39: computes `[~,eigvals]` using `[~,eigvals]=eig(M,'vector')`.
- Lines 42: computes `eigvals` using `eigvals=sort(eigvals,'ascend')`.
- Lines 45: computes `ax` using `ax=2*eigvals(3)-(eigvals(1)+eigvals(2))`.
- Lines 46: computes `rh` using `rh=eigvals(2)-eigvals(1)`.
- Lines 47: computes `iso` using `iso=mean(eigvals)`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(M)`. People say it would be terrible if we made all girls pretty. I think it would be great.
  - Representative operation: `if (~isnumeric(M))||(~isreal(M))|| (~all(size(M)==[3 3]))||(~issymmetric(M))`.
  - Representative operation: `(~all(size(M)==[3 3]))||(~issymmetric(M))`.

## Parameters / inputs

- M -a real symmetric 3x3 matrix

## Outputs

- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- ax -interaction axiality, defined as 2*zz-(xx+yy)
- in terms of eigenvalues
- rh -interaction rhombicity, defined as (yy-xx) in
- terms of eigenvalues
- eigvals -interaction tensor eigenvalues in Mehring order
- Note: eigenvalues [xx yy zz] are sorted in Mehring order, that
- is xx<=yy<=zz
- Note: Euler angles are not returned because the transformation
- in question is ill-defined

## Implementation structure

- Computes axiality and rhombicity of a symmetric 3x3 interaction
- tensor from the corresponding matrix. Syntax:
- [iso,ax,rh]=mat2axrh(M)
- M -a real symmetric 3x3 matrix
- iso -isotropic part of the interaction, defined as
- (xx+yy+zz)/3 in terms of eigenvaues
- ax -interaction axiality, defined as 2*zz-(xx+yy)
- in terms of eigenvalues
- rh -interaction rhombicity, defined as (yy-xx) in
- terms of eigenvalues
- eigvals -interaction tensor eigenvalues in Mehring order
- Note: eigenvalues [xx yy zz] are sorted in Mehring order, that

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `eigvals()`, `all()`, `issymmetric()`.
