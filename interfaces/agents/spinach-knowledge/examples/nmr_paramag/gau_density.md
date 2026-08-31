# examples/nmr_paramag/gau_density.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/gau_density.m`
- Signature: `gau_density()`
- Total lines: 30

## Purpose

Simple calculation of the PCS field of a Gaussian distribution of the electron probability density.

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Set problem dimensions; implemented by `dim=100; ext=10; sigma=0.5`.
- Lines 12-13: Get a 3D grid; implemented by `[X,Y,Z]=meshgrid(linspace(-ext,ext,dim),linspace(-ext,ext,dim),linspace(-ext,ext,dim))`.
- Lines 15-16: Pick a reasonable susceptibility tensor; implemented by `R=euler2dcm(pi/4,pi/5,pi/6)`.
- Lines 19-20: Get electron distribution; implemented by `probden=(1/sqrt((2*pi)^3*sigma^3))*exp(-(X.^2+Y.^2+Z.^2)/(2*sigma))`.
- Lines 22-23: Solve Kuprov equation; implemented by `[~,pcs_3d]=kpcs(probden,chi,[-ext ext -ext ext -ext ext],[0 0 0],'fft')`.
- Lines 25-26: Plot the solution; implemented by `kfigure(); volplot(pcs_3d,[-ext ext -ext ext -ext ext])`.

### Key state/data transformations

- Lines 10: computes `dim` using `dim=100; ext=10; sigma=0.5`.
- Lines 13: computes `[X,Y,Z]` using `[X,Y,Z]=meshgrid(linspace(-ext,ext,dim),linspace(-ext,ext,dim),linspace(-ext,ext,dim))`.
- Lines 16: computes `R` using `R=euler2dcm(pi/4,pi/5,pi/6)`.
- Lines 17: computes `chi` using `chi=R*[-0.1 0 0; 0 -0.2 0; 0 0 0.3]*R'`.
- Lines 20: computes `probden` using `probden=(1/sqrt((2*pi)^3*sigma^3))*exp(-(X.^2+Y.^2+Z.^2)/(2*sigma))`.
- Lines 23: computes `[~,pcs_3d]` using `[~,pcs_3d]=kpcs(probden,chi,[-ext ext -ext ext -ext ext],[0 0 0],'fft')`.

## Implementation structure

- Simple calculation of the PCS field of a Gaussian distribution of the
- electron probability density.
- Set problem dimensions
- Get a 3D grid
- Pick a reasonable susceptibility tensor
- Get electron distribution
- Solve Kuprov equation
- Plot the solution

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `kpcs()`, `kfigure()`, `volplot()`, `ktitle()`.
