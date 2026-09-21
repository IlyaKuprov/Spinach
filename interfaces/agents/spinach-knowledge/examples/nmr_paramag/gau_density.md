# examples/nmr_paramag/gau_density.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_paramag/gau_density.m`
- Signature: `gau_density()`
- Total lines: 30

## Purpose

Simple calculation of the PCS field of a Gaussian distribution of the electron probability density.

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.

## Numerical / algorithmic content

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
