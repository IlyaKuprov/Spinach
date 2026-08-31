# experiments/pseudocon/kpcs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/kpcs.m`
- Signature: `[pcs_vals,pcs_cube]=kpcs(probden,chi,ranges,nxyz,method)`
- Total lines: 122

## Purpose

Computes the three-dimensional distribution of pseudocontact shift field by solving Kuprov equation for PCS. Syntax: [pcs_cube,pcs_vals]=kpcs(probden,chi,ranges,nxyz)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(probden,chi,ranges,nxyz)`.
- Lines 41-42: Isolate rank 2 component of chi; implemented by `[~,~,rank2]=mat2sphten(chi)`.
- Lines 45-46: Compute axis extents; implemented by `extents(1)=ranges(2)-ranges(1)`.
- Lines 50-51: Decide the method; implemented by `switch method`.
- Lines 55-58: Generate Fourier derivative multipliers; implemented by `[X,Y,Z]=ndgrid(fftdiff(1,size(probden,1),extents(1)/size(probden,1)), fftdiff(1,size(probden,2),extents(2)/size(probden,2)), fftdiff(1,size(probden,3),extents(3)/size(pr…`.
- Lines 60-61: Get Laplace operator in Fourier space; implemented by `L=X.^2+Y.^2+Z.^2; L(L==0)=1`.
- Lines 63-66: Get Kuprov operator in Fourier space; implemented by `K=-(1/3)*(chi(1,1)*X.*X+chi(1,2)*X.*Y+chi(1,3)*X.*Z+ chi(2,1)*Y.*X+chi(2,2)*Y.*Y+chi(2,3)*Y.*Z+ chi(3,1)*Z.*X+chi(3,2)*Z.*Y+chi(3,3)*Z.*Z)`.
- Lines 68-69: Solve Kuprov equation (in ppm); implemented by `pcs_cube=1e6*real(ifftn(K.*fftn(probden)./L))`.
- Lines 73-74: Get 3-point Laplace operator; implemented by `L=fdlap(size(probden),extents,3)`.
- Lines 76-77: Get 3-point Kuprov operator; implemented by `K=fdkup(size(probden),extents,chi,3)`.
- Lines 79-80: Solve Kuprov equation (in ppm); implemented by `pcs_cube=1e6*reshape(cgs(L,K*probden(:),1e-10,1000),size(probden))`.
- Lines 84-85: Complain and bomb out; implemented by `error('unknown solver.')`.
- Lines 89-90: Deallocate variables; implemented by `clear('L','K')`.
- Lines 92-95: Compute PCS values at the nuclei; implemented by `[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(pcs_cube,1)), linspace(ranges(3),ranges(4),size(pcs_cube,2)), linspace(ranges(5),ranges(6),size(pcs_cube,3)))`.

### Control flow inferred from the code

- Line 51: dispatches on `method`; cases `'fft'`, `'fdiff'`.

### Key state/data transformations

- Lines 42: computes `[~,~,rank2]` using `[~,~,rank2]=mat2sphten(chi)`.
- Lines 43: computes `chi` using `chi=sphten2mat([],[],rank2)`.
- Lines 46: computes `extents(1)` using `extents(1)=ranges(2)-ranges(1)`.
- Lines 47: computes `extents(2)` using `extents(2)=ranges(4)-ranges(3)`.
- Lines 48: computes `extents(3)` using `extents(3)=ranges(6)-ranges(5)`.
- Lines 56-58: computes `[X,Y,Z]` using `[X,Y,Z]=ndgrid(fftdiff(1,size(probden,1),extents(1)/size(probden,1)), fftdiff(1,size(probden,2),extents(2)/size(probden,2)), fftdiff(1,size(probden,3),extents(3)/size(pr…`.
- Lines 64-66: computes `K` using `K=-(1/3)*(chi(1,1)*X.*X+chi(1,2)*X.*Y+chi(1,3)*X.*Z+ chi(2,1)*Y.*X+chi(2,2)*Y.*Y+chi(2,3)*Y.*Z+ chi(3,1)*Z.*X+chi(3,2)*Z.*Y+chi(3,3)*Z.*Z)`.
- Lines 69: computes `pcs_cube` using `pcs_cube=1e6*real(ifftn(K.*fftn(probden)./L))`.
- Lines 74: computes `L` using `L=fdlap(size(probden),extents,3)`.
- Lines 96: computes `pcs_vals` using `pcs_vals=interpn(X,Y,Z,pcs_cube,nxyz(:,1),nxyz(:,2),nxyz(:,3),'cubic')`.

### Local helper functions

- Line 101: `grumble()` — `function grumble(probden,chi,ranges,nxyz)`.
  - Representative operation: `if (~isnumeric(nxyz))||(size(nxyz,2)~=3)||(~isreal(nxyz))`.
  - Representative operation: `error('nxyz parameter should be a real matrix with three columns.')`.

## Parameters / inputs

- probden -unpaired electron probability density (not spin density) cube
- chi -electron magnetic susceptibility tensor in cubic Angstroms
- ranges -Cartesian axis extents for the unpaired electron probability
- density cube as [xmin xmax ymin ymax zmin zmax] in Angstroms
- nxyz -nuclear coordinates as [x y z] with multiple rows) at which
- PCS is to be evaluated, in Angstroms
- Output:
- pcs_vals -pseudocontact shift in ppm at each nucleus
- pcs_cube -pseudocontact shift field on the same grid as the unpaired
- electron probability density supplied
- Note: minimal three-point schemes are used for the finite difference
- operators. Increase stencil size if you have enough memory.
- Note: for further information on the equations and algorithms used in this
- function see http://dx.doi.org/10.1039/C4CP03106G

## Implementation structure

- Computes the three-dimensional distribution of pseudocontact shift field
- by solving Kuprov equation for PCS. Syntax:
- [pcs_cube,pcs_vals]=kpcs(probden,chi,ranges,nxyz)
- probden -unpaired electron probability density (not spin density) cube
- chi -electron magnetic susceptibility tensor in cubic Angstroms
- ranges -Cartesian axis extents for the unpaired electron probability
- density cube as [xmin xmax ymin ymax zmin zmax] in Angstroms
- nxyz -nuclear coordinates as [x y z] with multiple rows) at which
- PCS is to be evaluated, in Angstroms
- Output:
- pcs_vals -pseudocontact shift in ppm at each nucleus
- pcs_cube -pseudocontact shift field on the same grid as the unpaired

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mat2sphten()`, `sphten2mat()`, `extents()`, `ranges()`, `fftdiff()`, `chi()`, `fdlap()`, `fdkup()`, `cgs()`, `probden()`, `clear()`, `interpn()`, `nxyz()`, `ndims()`.
