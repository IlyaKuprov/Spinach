# examples/nmr_paramag/gau_density.m

- Signature: `gau_density()`

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
