# experiments/pseudocon/lpcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/lpcs.m`
- Signature: `theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)`
- Total lines: 138

## Purpose

Computes PCS from the multipole moments of the paramagnetic centre probability density as described in Equation 33 is used under the assumption that all nuclei are outside the bounding sphere shown in Figure 1 of the paper.

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 51-52: Check consistency; implemented by `grumble(nxyz,mxyz,ranks,Ilm,chi)`.
- Lines 54-55: Put the metal at the origin; implemented by `nxyz=nxyz-kron(mxyz,ones(size(nxyz,1),1))`.
- Lines 57-58: Form the susceptibility tensor; implemented by `if numel(chi)==5`.
- Lines 64-65: Compute the irreducible components; implemented by `[~,~,chi2]=qform2sph(chi)`.
- Lines 67-68: Build look-up tables; implemented by `A=zeros(max(ranks)+1,1)`.
- Lines 81-82: Convert coordinates to spherical; implemented by `[r,theta,phi]=xyz2sph(nxyz(:,1),nxyz(:,2),nxyz(:,3))`.
- Lines 84-85: Preallocate the answer; implemented by `theo_pcs=zeros(numel(r),1)`.
- Lines 99-100: Clean up and convert to ppm; implemented by `theo_pcs=1e6*real(theo_pcs)`.

### Control flow inferred from the code

- Line 58: conditional branch on `numel(chi)==5`.
- Line 69: `for` loop over `n=0:max(ranks)`.
- Line 73: `for` loop over `n=0:max(ranks)`.
- Line 74: `for` loop over `p=-n:n`.
- Line 75: `for` loop over `q=-2:2`.
- Line 88: `for` loop over `i=1:numel(ranks)`.
- Line 90: `for` loop over `m=-l:l`.
- Line 91: `for` loop over `mp=-2:2`.

### Key state/data transformations

- Lines 55: computes `nxyz` using `nxyz=nxyz-kron(mxyz,ones(size(nxyz,1),1))`.
- Lines 59: computes `chi` using `chi=[chi(1) chi(2) chi(3)`.
- Lines 65: computes `[~,~,chi2]` using `[~,~,chi2]=qform2sph(chi)`.
- Lines 68: computes `A` using `A=zeros(max(ranks)+1,1)`.
- Lines 70: computes `A(n+1)` using `A(n+1)=cg_fast(n+2,0,n,0,2,0)`.
- Lines 72: computes `B` using `B=zeros(max(ranks)+1,2*(max(ranks)+1)+1,5)`.
- Lines 76: computes `B(n+1,p+n+1,q+3)` using `B(n+1,p+n+1,q+3)=cg_fast(n+2,p+q,n,p,2,q)`.
- Lines 82: computes `[r,theta,phi]` using `[r,theta,phi]=xyz2sph(nxyz(:,1),nxyz(:,2),nxyz(:,3))`.
- Lines 85: computes `theo_pcs` using `theo_pcs=zeros(numel(r),1)`.
- Lines 89: computes `l` using `l=ranks(i)`.

### Local helper functions

- Line 105: `grumble()` — `function grumble(nxyz,mxyz,L,Ilm,chi)`.
  - Representative operation: `if (~isnumeric(nxyz))||(~isreal(nxyz))||any(~isfinite(nxyz(:)))||(size(nxyz,2)~=3)`.
  - Representative operation: `error('nxyz must be an Nx3 array of finite real atomic coordinates.')`.

## Syntax

```matlab
theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)
```

## Parameters / inputs

- nxyz -nuclear coordinates as [x y z] with multiple
- rows, at which PCS is to be evaluated, in
- Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z],
- in Angstroms.
- ranks -array of multipole ranks supplied in Ilm
- Ilm -{[],[]} array of numbers corresponding to
- the integrals
- Int[rho(r,theta,phi)*Y_lm(theta,phi)*r^l*d^3r]
- for L=0 Ilm=N/2/sqrt(pi)
- for L=1 Ilm=[imag(I11) I10 real(I11)]
- for L=2 Ilm=[imag(I22) imag(I21) I20 ...
- real(I21) real(I22)]
- et cetera.
- chi -the 3x3 matrix or the five independent elements
- of the susceptibility tensor in cubic Angstroms
- Output:
- theo_pcs -predicted pseudocontact shift (in ppm) at
- each of the nuclei.

## Implementation structure

- Computes PCS from the multipole moments of the paramagnetic
- centre probability density as described in
- Equation 33 is used under the assumption that all nuclei are
- outside the bounding sphere shown in Figure 1 of the paper.
- theo_pcs=lpcs(nxyz,mxyz,ranks,Ilm,chi)
- nxyz -nuclear coordinates as [x y z] with multiple
- rows, at which PCS is to be evaluated, in
- Angstroms.
- mxyz -paramagnetic centre coordinates as [x y z],
- in Angstroms.
- ranks -array of multipole ranks supplied in Ilm
- Ilm -{[],[]} array of numbers corresponding to

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `chi()`, `qform2sph()`, `cg_fast()`, `xyz2sph()`, `nxyz()`, `ranks()`, `sign()`, `chi2()`, `spher_harmon()`, `any()`, `mxyz()`, `isvector()`, `isequal()`, `iscell()`.
