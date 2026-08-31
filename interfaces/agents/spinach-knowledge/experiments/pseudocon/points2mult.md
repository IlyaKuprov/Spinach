# experiments/pseudocon/points2mult.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/points2mult.m`
- Signature: `Ilm=points2mult(xyz,mxyz,rho,L,method)`
- Total lines: 107

## Purpose

Computes multipole moments from a set of points with user-specified spin populations. Syntax: Ilm=points2mult(xyz,mxyz,rho,L,method) The multipoles in question are described in Equation 32 of:

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(xyz,mxyz,rho,L)`.
- Lines 43-44: Put the metal at the origin; implemented by `xyz=xyz-kron(mxyz,ones(size(xyz,1),1))`.
- Lines 46-47: Convert coordinates to spherical; implemented by `[r,theta,phi]=xyz2sph(xyz(:,1),xyz(:,2),xyz(:,3))`.
- Lines 49-50: Compute the coefficients; implemented by `for n=1:numel(L)`.
- Lines 62-63: Process the options; implemented by `switch method`.
- Lines 67-68: Compute the volume element; implemented by `x=unique(xyz(:,1)); dx=(max(x)-min(x))/(numel(x)-1)`.
- Lines 72-73: Apply the volume element; implemented by `Ilm=Ilm*dx*dy*dz`.
- Lines 77-79: No action necessary; implemented by `end`.

### Control flow inferred from the code

- Line 50: `for` loop over `n=1:numel(L)`.
- Line 52: `for` loop over `m=0:l`.
- Line 53: conditional branch on `m==0`.
- Line 63: dispatches on `method`; cases `'grid'`, `'points'`.

### Key state/data transformations

- Lines 44: computes `xyz` using `xyz=xyz-kron(mxyz,ones(size(xyz,1),1))`.
- Lines 47: computes `[r,theta,phi]` using `[r,theta,phi]=xyz2sph(xyz(:,1),xyz(:,2),xyz(:,3))`.
- Lines 51: computes `l` using `l=L(n)`.
- Lines 54: computes `Ilm{n}(l+1)` using `Ilm{n}(l+1)=sum(rho.*r.^l.*spher_harmon(l,0,theta,phi))`.
- Lines 56: computes `Ilm{n}( m+l+1)` using `Ilm{n}( m+l+1)=+real(sum(rho.*r.^l.*spher_harmon(l,m,theta,phi)))`.
- Lines 57: computes `Ilm{n}(-m+l+1)` using `Ilm{n}(-m+l+1)=-imag(sum(rho.*r.^l.*spher_harmon(l,m,theta,phi)))`.
- Lines 68: computes `x` using `x=unique(xyz(:,1)); dx=(max(x)-min(x))/(numel(x)-1)`.
- Lines 69: computes `y` using `y=unique(xyz(:,2)); dy=(max(y)-min(y))/(numel(y)-1)`.
- Lines 70: computes `z` using `z=unique(xyz(:,3)); dz=(max(z)-min(z))/(numel(z)-1)`.
- Lines 73: computes `Ilm` using `Ilm=Ilm*dx*dy*dz`.

### Local helper functions

- Line 84: `grumble()` — `function grumble(xyz,mxyz,rho,L)`.
  - Representative operation: `if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,2)~=3)`.
  - Representative operation: `error('xyz must be an Nx3 array of atomic coordinates.')`.

## Parameters / inputs

- xyz -coordinates as [x y z] with multiple rows,
- at which density rho is evaluated, in Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z], in
- Angstroms.
- L -array of ranks of spherical harmonics of the
- probability density
- rho -column of the densities at the points xyz
- method -'points' or 'grid'. If the spin density is supplied
- as Mulliken spin populations at individual nuclei,
- choose 'points'; if the spin density is supplied as
- a probability on a uniform cubic grid obtained from
- ndgrid() function and vectorised, use 'grid'.
- Output:
- Ilm -multipole moments of the probability density

## Implementation structure

- Computes multipole moments from a set of points with user-specified
- spin populations. Syntax:
- Ilm=points2mult(xyz,mxyz,rho,L,method)
- The multipoles in question are described in Equation 32 of:
- xyz -coordinates as [x y z] with multiple rows,
- at which density rho is evaluated, in Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z], in
- Angstroms.
- L -array of ranks of spherical harmonics of the
- probability density
- rho -column of the densities at the points xyz
- method -'points' or 'grid'. If the spin density is supplied

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xyz2sph()`, `xyz()`, `spher_harmon()`.
