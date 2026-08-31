# experiments/pseudocon/ppcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/ppcs.m`
- Signature: `pcs=ppcs(nxyz,sxyz,chi)`
- Total lines: 80

## Purpose

Computes pseudocontact shift from a point electron centre at the nuclear coordinates supplied. Syntax: pred_pcs=ppcs(nxyz,mxyz,chi)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Form susceptibility tensor; implemented by `if numel(chi)==5`.
- Lines 39-40: Check consistency; implemented by `grumble(nxyz,sxyz,chi)`.
- Lines 42-43: Put the metal at the origin; implemented by `nxyz=nxyz-kron(sxyz,ones(size(nxyz,1),1))`.
- Lines 45-46: Convert coordinates to spherical; implemented by `[r,theta,phi]=xyz2sph(nxyz(:,1),nxyz(:,2),nxyz(:,3))`.
- Lines 48-49: Get the irreducible components; implemented by `[~,~,chi2]=qform2sph(chi)`.
- Lines 51-52: Preallocate the answer; implemented by `pcs=zeros(numel(r),1)`.
- Lines 54-55: Sum the spherical harmonic expansion; implemented by `for m=[2 1 0 -1 -2]`.
- Lines 59-60: Convert to ppm and clean up; implemented by `pcs=1e6*real(pcs)`.

### Control flow inferred from the code

- Line 33: conditional branch on `numel(chi)==5`.
- Line 55: `for` loop over `m=[2 1 0 -1 -2]`.

### Key state/data transformations

- Lines 34: computes `chi` using `chi=[chi(1) chi(2) chi(3)`.
- Lines 43: computes `nxyz` using `nxyz=nxyz-kron(sxyz,ones(size(nxyz,1),1))`.
- Lines 46: computes `[r,theta,phi]` using `[r,theta,phi]=xyz2sph(nxyz(:,1),nxyz(:,2),nxyz(:,3))`.
- Lines 49: computes `[~,~,chi2]` using `[~,~,chi2]=qform2sph(chi)`.
- Lines 52: computes `pcs` using `pcs=zeros(numel(r),1)`.

### Local helper functions

- Line 65: `grumble()` — `function grumble(nxyz,sxyz,chi)`.
  - Representative operation: `if (~isnumeric(nxyz))||(size(nxyz,2)~=3)||(~isreal(nxyz))`.
  - Representative operation: `error('nxyz must be a real matrix with three columns.')`.

## Parameters / inputs

- chi -magnetic susceptibility tensor in cubic Angstroms
- as a 3x3 matrix, or its five unique components
- ordered as
- [chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is to be evaluated, in Angstroms.
- sxyz -susceptibility centre coordinates as [x y z], in
- Angstroms.
- Output:
- pcs -predicted pseudocontact shift (in ppm) at each of
- the nuclei.

## Implementation structure

- Computes pseudocontact shift from a point electron centre at the
- nuclear coordinates supplied. Syntax:
- pred_pcs=ppcs(nxyz,mxyz,chi)
- chi -magnetic susceptibility tensor in cubic Angstroms
- as a 3x3 matrix, or its five unique components
- ordered as
- [chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]
- nxyz -nuclear coordinates as [x y z] with multiple rows,
- at which PCS is to be evaluated, in Angstroms.
- sxyz -susceptibility centre coordinates as [x y z], in
- Angstroms.
- Output:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `chi()`, `grumble()`, `xyz2sph()`, `nxyz()`, `qform2sph()`, `chi2()`, `spher_harmon()`.
